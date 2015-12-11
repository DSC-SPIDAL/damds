package edu.indiana.soic.spidal.damds;

import com.google.common.base.Optional;
import com.google.common.base.Stopwatch;
import com.google.common.base.Strings;
import edu.indiana.soic.spidal.common.*;
import edu.indiana.soic.spidal.configuration.ConfigurationMgr;
import edu.indiana.soic.spidal.configuration.section.DAMDSSection;
import edu.indiana.soic.spidal.damds.timing.*;
import mpi.MPI;
import mpi.MPIException;
import net.openhft.lang.io.Bytes;
import org.apache.commons.cli.*;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.LongBuffer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;
import java.util.stream.IntStream;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;


public class Program {
    private static Options programOptions = new Options();

    static {
        programOptions.addOption(
            String.valueOf(Constants.CMD_OPTION_SHORT_C),
            Constants.CMD_OPTION_LONG_C, true,
            Constants.CMD_OPTION_DESCRIPTION_C);
        programOptions.addOption(
            String.valueOf(Constants.CMD_OPTION_SHORT_N),
            Constants.CMD_OPTION_LONG_N, true,
            Constants.CMD_OPTION_DESCRIPTION_N);
        programOptions.addOption(
            String.valueOf(Constants.CMD_OPTION_SHORT_T),
            Constants.CMD_OPTION_LONG_T, true,
            Constants.CMD_OPTION_DESCRIPTION_T);

        programOptions.addOption(Constants.CMD_OPTION_SHORT_MMAPS, true, Constants.CMD_OPTION_DESCRIPTION_MMAPS);
        programOptions.addOption(
            Constants.CMD_OPTION_SHORT_MMAP_SCRATCH_DIR, true,
            Constants.CMD_OPTION_DESCRIPTION_MMAP_SCRATCH_DIR);
    }

    // Constants
    public static final double INV_SHORT_MAX = 1.0 / Short.MAX_VALUE;
    public static final double SHORT_MAX = Short.MAX_VALUE;

    // Calculated Constants
    public static double INV_SUM_OF_SQUARE;

    // Arrays
    public static double[][] preX;
    public static double[][] BC;
    public static double[][] MMr;
    public static double[][] MMAp;

    public static double[][][] threadPartialBofZ;
    public static double[][][] threadPartialMM;

    public static double[] partialSigma;
    public static double[][] vArray;

    //Config Settings
    public static DAMDSSection config;
    public static ByteOrder byteOrder;
    public static short[][] distances;
    public static WeightsWrap weights;

    public static int BlockSize;

    /**
     * Weighted SMACOF based on Deterministic Annealing algorithm
     *
     * @param args command line arguments to the program, which should include
     *             -c path to config file
     *             -t number of threads
     *             -n number of nodes
     *             The options may also be given as longer names
     *             --configFile, --threadCount, and --nodeCount respectively
     */
    public static void  main(String[] args) {
        Optional<CommandLine> parserResult =
            parseCommandLineArguments(args, programOptions);

        if (!parserResult.isPresent()) {
            System.out.println(Constants.ERR_PROGRAM_ARGUMENTS_PARSING_FAILED);
            new HelpFormatter()
                .printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        CommandLine cmd = parserResult.get();
        if (!(cmd.hasOption(Constants.CMD_OPTION_LONG_C) &&
              cmd.hasOption(Constants.CMD_OPTION_LONG_N) &&
              cmd.hasOption(Constants.CMD_OPTION_LONG_T))) {
            System.out.println(Constants.ERR_INVALID_PROGRAM_ARGUMENTS);
            new HelpFormatter()
                .printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        try {
            //  Read Metadata using this as source of other metadata
            readConfiguration(cmd);

            //  Set up MPI and threads parallelism
            ParallelOps.setupParallelism(args);
            // Note - a barrier to get cleaner timings
            ParallelOps.worldProcsComm.barrier();
            Stopwatch mainTimer = Stopwatch.createStarted();

            ParallelOps.setParallelDecomposition(
                config.numberDataPoints, config.targetDimension);
            initializeTimers();

            Utils.printMessage("\n== DAMDS run started on " + new Date() + " ==\n");
            Utils.printMessage(config.toString(false));

            readDistancesAndWeights(config.isSammon);
            RefObj<Integer> missingDistCount = new RefObj<>();
            DoubleStatistics distanceSummary = calculateStatistics(
                distances, weights, missingDistCount);
            double missingDistPercent = missingDistCount.getValue() /
                                        (Math.pow(config.numberDataPoints, 2));
            INV_SUM_OF_SQUARE = 1.0/distanceSummary.getSumOfSquare();
            Utils.printMessage(
                "\nDistance summary... \n" + distanceSummary.toString() +
                "\n  MissingDistPercentage=" +
                missingDistPercent);

            weights.setAvgDistForSammon(distanceSummary.getAverage());
            changeZeroDistancesToPostiveMin(distances, distanceSummary.getPositiveMin());

            // Allocating point arrays once for all
            allocateArrays();

            if (Strings.isNullOrEmpty(config.initialPointsFile)){
                generateInitMapping(
                    config.numberDataPoints, config.targetDimension, preX);
            } else {
                readInitMapping(config.initialPointsFile, preX);
            }
            double tCur = 0.0;
            double tMax = distanceSummary.getMax() / Math.sqrt(2.0 * config.targetDimension);
            double tMin = config.tMinFactor * distanceSummary.getPositiveMin() / Math.sqrt(2.0 * config.targetDimension);

            generateVArray(distances, weights, vArray);
            double preStress = calculateStress(
                preX, config.targetDimension, tCur, distances, weights,
                INV_SUM_OF_SQUARE, partialSigma);
            Utils.printMessage("\nInitial stress=" + preStress);

            tCur = config.alpha * tMax;

            mainTimer.stop();
            Utils.printMessage(
                "\nUp to the loop took " + mainTimer.elapsed(
                    TimeUnit.SECONDS) + " seconds");
            mainTimer.start();

            Stopwatch loopTimer = Stopwatch.createStarted();

            int loopNum = 0;
            double diffStress;
            double stress = -1.0;
            RefObj<Integer> outRealCGIterations = new RefObj<>(0);
            RefObj<Integer> cgCount = new RefObj<>(0);
            int smacofRealIterations = 0;
            while (true) {

                TemperatureLoopTimings.startTiming(
                    TemperatureLoopTimings.TimingTask.PRE_STRESS);
                preStress = calculateStress(
                    preX, config.targetDimension, tCur, distances, weights,
                    INV_SUM_OF_SQUARE, partialSigma);
                TemperatureLoopTimings.endTiming(
                    TemperatureLoopTimings.TimingTask.PRE_STRESS);

                diffStress = config.threshold + 1.0;

                Utils.printMessage(
                    String.format(
                        "\nStart of loop %d Temperature (T_Cur) %.5g",
                        loopNum, tCur));

                int itrNum = 0;
                cgCount.setValue(0);
                TemperatureLoopTimings.startTiming(
                    TemperatureLoopTimings.TimingTask.STRESS_LOOP);
                while (diffStress >= config.threshold) {

                    zeroOutArray(threadPartialMM);
                    StressLoopTimings.startTiming(
                        StressLoopTimings.TimingTask.BC);
                    calculateBC(
                        preX, config.targetDimension, tCur, distances,
                        weights, BlockSize, BC, threadPartialBofZ,
                        threadPartialMM);
                    StressLoopTimings.endTiming(
                        StressLoopTimings.TimingTask.BC);
                    // This barrier was necessary for correctness when using
                    // a single mmap file
                    ParallelOps.worldProcsComm.barrier();


                    StressLoopTimings.startTiming(
                        StressLoopTimings.TimingTask.CG);
                    calculateConjugateGradient(preX, config.targetDimension,
                                               config.numberDataPoints,
                                               BC,
                                               config.cgIter,
                                               config.cgErrorThreshold, cgCount,
                                               outRealCGIterations, weights,
                                               BlockSize, vArray, MMr, MMAp, threadPartialMM);
                    StressLoopTimings.endTiming(
                        StressLoopTimings.TimingTask.CG);


                    StressLoopTimings.startTiming(
                        StressLoopTimings.TimingTask.STRESS);
                    stress = calculateStress(
                        preX, config.targetDimension, tCur, distances, weights,
                        INV_SUM_OF_SQUARE, partialSigma);
                    StressLoopTimings.endTiming(
                        StressLoopTimings.TimingTask.STRESS);


                    diffStress = preStress - stress;
                    preStress = stress;

                    if ((itrNum % 10 == 0) || (itrNum >= config.stressIter)) {
                        Utils.printMessage(
                            String.format(
                                "  Loop %d Iteration %d Avg CG count %.5g " +
                                "Stress " +
                                "%.5g", loopNum, itrNum,
                                (cgCount.getValue() * 1.0 / (itrNum + 1)),
                                stress));
                    }
                    ++itrNum;
                    ++smacofRealIterations;
                }
                TemperatureLoopTimings.endTiming(
                    TemperatureLoopTimings.TimingTask.STRESS_LOOP);

                --itrNum;
                if (itrNum >=0 && !(itrNum % 10 == 0) && !(itrNum >=
                                                           config.stressIter)) {
                    Utils.printMessage(
                        String.format(
                            "  Loop %d Iteration %d Avg CG count %.5g " +
                            "Stress %.5g",
                            loopNum, itrNum,
                            (cgCount.getValue() * 1.0 / (itrNum + 1)), stress));
                }

                Utils.printMessage(
                    String.format(
                        "End of loop %d Total Iterations %d Avg CG count %.5g" +
                        " Stress %.5g",
                        loopNum, (itrNum + 1),
                        (cgCount.getValue() * 1.0 / (itrNum + 1)), stress));

                if (tCur == 0)
                    break;
                tCur *= config.alpha;
                if (tCur < tMin)
                    tCur = 0;
                ++loopNum;
            }
            loopTimer.stop();

            double QoR1 = stress / (config.numberDataPoints * (config.numberDataPoints - 1) / 2);
            double QoR2 = QoR1 / (distanceSummary.getAverage() * distanceSummary.getAverage());

            Utils.printMessage(
                String.format(
                    "Normalize1 = %.5g Normalize2 = %.5g", QoR1, QoR2));
            Utils.printMessage(
                String.format(
                    "Average of Delta(original distance) = %.5g",
                    distanceSummary.getAverage()));


             // TODO Fix error handling here
            /*if (Strings.isNullOrEmpty(config.labelFile) || config.labelFile.toUpperCase().endsWith(
                "NOLABEL")) {
                try {
                    writeOuput(preX, config.pointsFile);
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
            } else {
                try {
                    writeOuput(preX, config.labelFile, config.pointsFile);
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
            }*/

            Double finalStress = calculateStress(
                preX, config.targetDimension, tCur, distances, weights,
                INV_SUM_OF_SQUARE, partialSigma);

            mainTimer.stop();


            Utils.printMessage("Finishing DAMDS run ...");
            long totalTime = mainTimer.elapsed(TimeUnit.MILLISECONDS);
            long temperatureLoopTime = loopTimer.elapsed(TimeUnit.MILLISECONDS);
            Utils.printMessage(
                String.format(
                    "  Total Time: %s (%d ms) Loop Time: %s (%d ms)",
                    formatElapsedMillis(totalTime), totalTime,
                    formatElapsedMillis(temperatureLoopTime), temperatureLoopTime));
            Utils.printMessage("  Total Loops: " + loopNum);
            Utils.printMessage("  Total Iterations: " + smacofRealIterations);
            Utils.printMessage(
                String.format(
                    "  Total CG Iterations: %d Avg. CG Iterations: %.5g",
                    outRealCGIterations.getValue(), (outRealCGIterations.getValue() * 1.0) / smacofRealIterations));
            Utils.printMessage("  Final Stress:\t" + finalStress);

            printTimings(totalTime, temperatureLoopTime);

            Utils.printMessage("== DAMDS run completed on " + new Date() + " ==");

            ParallelOps.tearDownParallelism();
        }
        catch (MPIException | IOException | InterruptedException e) {
            Utils.printAndThrowRuntimeException(new RuntimeException(e));
        }
    }

    private static void allocateArrays() {
        // Allocating point arrays once for all
        final int numberDataPoints = config.numberDataPoints;
        final int targetDimension = config.targetDimension;
        preX = new double[numberDataPoints][targetDimension];
        BC = new double[numberDataPoints][targetDimension];
        MMr = new double[numberDataPoints][targetDimension];
        MMAp = new double[numberDataPoints][targetDimension];
        final int threadCount = ParallelOps.threadCount;
        threadPartialBofZ = new double[threadCount][][];
        threadPartialMM = new double[threadCount][][];
        vArray = new double[threadCount][];
        int threadRowCount;
        for (int i = 0; i < threadCount; ++i){
            threadRowCount = ParallelOps.threadRowCounts[i];
            threadPartialBofZ[i] = new double[threadRowCount][ParallelOps.globalColCount];
            threadPartialMM[i] = new double[threadRowCount][config.targetDimension];
            vArray[i] = new double[threadRowCount];
        }
        partialSigma = new double[threadCount];
    }

    private static void zeroOutArray(double[][] a){
        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> Arrays.fill(a[threadIdx], 0.0d)));
        }
        else {
            Arrays.fill(a[0],0.0d);
        }
    }

    private static void zeroOutArray(double[][][] a) {
        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> zeroOutArrayInternal(
                        ParallelOps.threadRowCounts[threadIdx],
                        a[threadIdx])));
        }
        else {
            zeroOutArrayInternal(ParallelOps.threadRowCounts[0], a[0]);
        }
    }

    private static void zeroOutArrayInternal(int threadRowCount, double[][] a){
        for (int j = 0; j < threadRowCount; ++j){
            Arrays.fill(a[j], 0.0d);
        }
    }


    private static void changeZeroDistancesToPostiveMin(
        short[][] distances, double positiveMin) {
        double tmpD;
        for (short[] distanceRow : distances){
            for (int j = 0; j < distanceRow.length; ++j){
                tmpD = distanceRow[j] * INV_SHORT_MAX;
                if (tmpD < positiveMin && tmpD >= 0.0){
                    distanceRow[j] = (short)(positiveMin * SHORT_MAX);
                }
            }
        }
    }

    private static void printTimings(long totalTime, long temperatureLoopTime)
        throws MPIException {
        String mainHeader =
            "  Total\tTempLoop\tPreStress\tStressLoop\tBC\tCG\tStress" +
            "\tBCInternal\tBComm\tBCInternalBofZ\tBCInternalMM\tCGMM" +
            "\tCGInnerProd\tCGLoop\tCGLoopMM\tCGLoopInnerProdPAP" +
            "\tCGLoopInnerProdR\tMMInternal\tMMComm\tBCMerge\tBCExtract\tMMMerge\tMMExtract\tStressInternal\tStressComm\tStressInternalComp";
        Utils.printMessage(
            mainHeader);
        String mainTimings = "  " + totalTime + '\t' + temperatureLoopTime +
                             '\t' +
                     TemperatureLoopTimings.getAverageTime(
                         TemperatureLoopTimings.TimingTask.PRE_STRESS) + '\t' +
                     TemperatureLoopTimings.getAverageTime(
                         TemperatureLoopTimings.TimingTask.STRESS_LOOP) + '\t' +
                     StressLoopTimings.getAverageTime(
                         StressLoopTimings.TimingTask.BC) + '\t' +
                     StressLoopTimings.getAverageTime(
                         StressLoopTimings.TimingTask.CG) + '\t' +
                     StressLoopTimings.getAverageTime(
                         StressLoopTimings.TimingTask.STRESS) + '\t' +
                     BCTimings.getAverageTime(
                         BCTimings.TimingTask.BC_INTERNAL) + '\t' +
                     BCTimings.getAverageTime(
                         BCTimings.TimingTask.COMM) + '\t' +
                     BCInternalTimings.getAverageTime(
                         BCInternalTimings.TimingTask.BOFZ) + '\t' +
                     BCInternalTimings.getAverageTime(
                         BCInternalTimings.TimingTask.MM) + '\t' +
                     CGTimings.getAverageTime(
                         CGTimings.TimingTask.MM) + '\t' +
                     CGTimings.getAverageTime(
                         CGTimings.TimingTask.INNER_PROD) + '\t' +
                     CGTimings.getAverageTime(
                         CGTimings.TimingTask.CG_LOOP) + '\t' +
                     CGLoopTimings.getAverageTime(
                         CGLoopTimings.TimingTask.MM) + '\t' +
                     CGLoopTimings.getAverageTime(
                         CGLoopTimings.TimingTask.INNER_PROD_PAP) + '\t' +
                     CGLoopTimings.getAverageTime(
                         CGLoopTimings.TimingTask.INNER_PROD_R) + '\t' +
                     MMTimings.getAverageTime(
                         MMTimings.TimingTask.MM_INTERNAL) + '\t' +
                     MMTimings.getAverageTime(
                         MMTimings.TimingTask.COMM) + '\t' +
                     BCTimings.getAverageTime(BCTimings.TimingTask.BC_MERGE) + '\t' +
                     BCTimings.getAverageTime(BCTimings.TimingTask.BC_EXTRACT) + '\t' +
                     MMTimings.getAverageTime(MMTimings.TimingTask.MM_MERGE) + '\t' +
                     MMTimings.getAverageTime(MMTimings.TimingTask.MM_EXTRACT) + '\t' +
                     StressTimings.getAverageTime(
                         StressTimings.TimingTask.STRESS_INTERNAL) + '\t' +
                     StressTimings.getAverageTime(
                           StressTimings.TimingTask.COMM) + '\t' +
                     StressInternalTimings.getAverageTime(
                           StressInternalTimings.TimingTask.COMP);

        Utils.printMessage(
            mainTimings);

        // Accumulated times as percentages of the total time
        // Note, for cases where threads are present, the max time
        // is taken as the total time of that component for the particular rank
        String percentHeader =
            "  Total%\tTempLoop%\tPreStress%\tStressLoop%\tBC\tCG%\tStress%" +
            "\tBCInternal%\tBComm%\tBCInternalBofZ%\tBCInternalMM%\tCGMM%" +
            "\tCGInnerProd%\tCGLoop%\tCGLoopMM%\tCGLoopInnerProdPAP%" +
            "\tCGLoopInnerProdR%\tMMInternal%\tMMComm%\tBCMerge%\tBCExtract%"
            + "\tMMMerge%\tMMExtract%";
        Utils.printMessage(
            percentHeader);
        String percentTimings =
            "  " + 1.0 + '\t' + (temperatureLoopTime *1.0 / totalTime) + '\t' +
            TemperatureLoopTimings.getTotalTime(
                TemperatureLoopTimings.TimingTask.PRE_STRESS) * 1.0 / totalTime +
            '\t' +
            TemperatureLoopTimings.getTotalTime(
                TemperatureLoopTimings.TimingTask.STRESS_LOOP) * 1.0 / totalTime +
            '\t' +
            StressLoopTimings.getTotalTime(
                StressLoopTimings.TimingTask.BC) * 1.0 / totalTime + '\t' +
            StressLoopTimings.getTotalTime(
                StressLoopTimings.TimingTask.CG) * 1.0 / totalTime + '\t' +
            StressLoopTimings.getTotalTime(
                StressLoopTimings.TimingTask.STRESS) * 1.0 / totalTime + '\t' +
            BCTimings.getTotalTime(
                BCTimings.TimingTask.BC_INTERNAL) * 1.0 / totalTime + '\t' +
            BCTimings.getTotalTime(
                BCTimings.TimingTask.COMM) * 1.0 / totalTime + '\t' +
            BCInternalTimings.getTotalTime(
                BCInternalTimings.TimingTask.BOFZ) * 1.0 / totalTime + '\t' +
            BCInternalTimings.getTotalTime(
                BCInternalTimings.TimingTask.MM) * 1.0 / totalTime + '\t' +
            CGTimings.getTotalTime(
                CGTimings.TimingTask.MM) * 1.0 / totalTime + '\t' +
            CGTimings.getTotalTime(
                CGTimings.TimingTask.INNER_PROD) * 1.0 / totalTime + '\t' +
            CGTimings.getTotalTime(
                CGTimings.TimingTask.CG_LOOP) * 1.0 / totalTime + '\t' +
            CGLoopTimings.getTotalTime(
                CGLoopTimings.TimingTask.MM) * 1.0 / totalTime + '\t' +
            CGLoopTimings.getTotalTime(
                CGLoopTimings.TimingTask.INNER_PROD_PAP) * 1.0 / totalTime + '\t' +
            CGLoopTimings.getTotalTime(
                CGLoopTimings.TimingTask.INNER_PROD_R) * 1.0 / totalTime + '\t' +
            MMTimings.getTotalTime(
                MMTimings.TimingTask.MM_INTERNAL) * 1.0 / totalTime + '\t' +
            MMTimings.getTotalTime(
                MMTimings.TimingTask.COMM) * 1.0 / totalTime + '\t' +
            BCTimings.getTotalTime(BCTimings.TimingTask.BC_MERGE) * 1.0 / totalTime + '\t' +
            BCTimings.getTotalTime(BCTimings.TimingTask.BC_EXTRACT) * 1.0 / totalTime + '\t' +
            MMTimings.getTotalTime(MMTimings.TimingTask.MM_MERGE) * 1.0 / totalTime + '\t' +
            MMTimings.getTotalTime(MMTimings.TimingTask.MM_EXTRACT) * 1.0 / totalTime + '\t' +
            StressTimings.getTotalTime(
                StressTimings.TimingTask.STRESS_INTERNAL) * 1.0 / totalTime + '\t' +
            StressTimings.getTotalTime(
                StressTimings.TimingTask.COMM) * 1.0 / totalTime + '\t' +
            StressInternalTimings.getTotalTime(
                StressInternalTimings.TimingTask.COMP) * 1.0 / totalTime;
        Utils.printMessage(
            percentTimings);

        // Timing (total timings) distribution against rank/threadID for,
        // 1. TempLoop time (MPI only)
        // 2. Stress (within TempLoop) (MPI only)
        // 3. BCComm (MPI only)
        // 4. MMComm (MPI only)
        // 5. BCInternalBofZ (has MPI+threads)
        // 6. BCInternalMM (has MPI+threads)
        // 7. MMInternal (has MPI+threads)

        long[] temperatureLoopTimeDistribution =
            getTemperatureLoopTimeDistribution(temperatureLoopTime);
        long[] stressTimeDistribution = StressLoopTimings
            .getTotalTimeDistribution(StressLoopTimings.TimingTask.STRESS);
        long[] bcCommTimeDistribution = BCTimings.getTotalTimeDistribution(
            BCTimings.TimingTask.COMM);
        long[] bcInternalBofZTimeDistribution =
            BCInternalTimings.getTotalTimeDistribution(
                BCInternalTimings.TimingTask.BOFZ);
        long[] bcInternalMMTimeDistribution =
            BCInternalTimings.getTotalTimeDistribution(
                BCInternalTimings.TimingTask.MM);
        long[] mmInternalTimeDistribution = MMTimings.getTotalTimeDistribution(
            MMTimings.TimingTask.MM_INTERNAL);
        long[] mmCommTimeDistribution = MMTimings.getTotalTimeDistribution(
            MMTimings.TimingTask.COMM);
        long[] bcMergeTimeDistribution = BCTimings.getTotalTimeDistribution(
            BCTimings.TimingTask.BC_MERGE);
        long[] bcExtractTimeDistribution = BCTimings.getTotalTimeDistribution(
            BCTimings.TimingTask.BC_EXTRACT);
        long[] mmMergeTimeDistribution = MMTimings.getTotalTimeDistribution(
            MMTimings.TimingTask.MM_MERGE);
        long[] mmExtractTimeDistribution = MMTimings.getTotalTimeDistribution(
            MMTimings.TimingTask.MM_EXTRACT);
        long[] stressCommTimeDistribution = StressTimings.getTotalTimeDistribution(
            StressTimings.TimingTask.COMM);
        long[] stressInternalCompTimeDistribution =
            StressInternalTimings.getTotalTimeDistribution(
                StressInternalTimings.TimingTask.COMP);

        // Count distributions
        long[] temperatureLoopCountDistribution =
            TemperatureLoopTimings.getCountDistribution(
                TemperatureLoopTimings.TimingTask.PRE_STRESS);
        long[] stressLoopCountDistribution =
            StressLoopTimings.getCountDistribution(
                StressLoopTimings.TimingTask.BC);
        long[] cgLoopCountDistribution = CGLoopTimings.getCountDistribution(
            CGLoopTimings.TimingTask.INNER_PROD_PAP);
        long[] mmInternalCountDistribution = MMTimings.getCountDistribution(
            MMTimings.TimingTask.MM_INTERNAL);

        if (ParallelOps.worldProcRank == 0) {
            try (BufferedWriter writer = Files.newBufferedWriter(
                Paths.get(config.timingFile), StandardOpenOption.WRITE,
                StandardOpenOption.CREATE)) {
                PrintWriter printWriter = new PrintWriter(writer, true);
                printWriter.println(mainHeader.trim());
                printWriter.println(mainTimings.trim());
                printWriter.println();
                printWriter.println(percentHeader.trim());
                printWriter.println(percentTimings.trim());
                printWriter.println();

                prettyPrintArray("Temperature Loop Timing Distribution",
                                 temperatureLoopTimeDistribution, printWriter);
                prettyPrintArray("Stress Timing Distribution", stressTimeDistribution, printWriter);
                prettyPrintArray("BCComm Timing Distribution",
                                 bcCommTimeDistribution, printWriter);
                prettyPrintArray("MMComm Timing Distribution",
                                 mmCommTimeDistribution, printWriter);
                prettyPrintArray("StressComm Timing Distribution",
                    stressCommTimeDistribution, printWriter);
                prettyPrintArray("BCInternalBofZ Timing Distribution",
                                 bcInternalBofZTimeDistribution, printWriter);
                prettyPrintArray("BCInternalMM Timing Distribution",
                                 bcInternalMMTimeDistribution, printWriter);
                prettyPrintArray("MMInternal Timing Distribution",
                                 mmInternalTimeDistribution, printWriter);
                prettyPrintArray("StressInternalComp Timing Distribution",
                    stressInternalCompTimeDistribution, printWriter);
                prettyPrintArray("BC Merge Timing Distribution",
                                 bcMergeTimeDistribution, printWriter);
                prettyPrintArray("BC Extract Timing Distribution",
                                 bcExtractTimeDistribution, printWriter);
                prettyPrintArray("MM Merge Timing Distribution", mmMergeTimeDistribution, printWriter);
                prettyPrintArray("MM Extract Timing Distribution", mmExtractTimeDistribution, printWriter);


                printWriter.println();

                prettyPrintArray("Temperature Loop Count Distribution",
                                 temperatureLoopCountDistribution, printWriter);
                prettyPrintArray("Stress Loop Count Distribution",
                                 stressLoopCountDistribution, printWriter);
                prettyPrintArray("CG Loop Count Distribution",
                                 cgLoopCountDistribution, printWriter);
                prettyPrintArray("MM Internal Count Distribution", mmInternalCountDistribution, printWriter);

                printWriter.println();
                // Print MPI rank
                String s = "";
                for (int i = 0; i < ParallelOps.worldProcsCount * ParallelOps.threadCount; ++i){
                    s += (i / ParallelOps.threadCount) + "\t";
                }
                printWriter.println(s);
                printWriter.println();

                printWriter.flush();
                printWriter.close();
            }
            catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    private static void prettyPrintArray(
        String title, long[] mmExtractTimeDistribution, PrintWriter printWriter) {
        String str;
        printWriter.println(title);
        str = Arrays.toString(mmExtractTimeDistribution);
        printWriter.println(
            str.substring(1, str.length() - 1).replace(',', '\t'));
        printWriter.println();
    }

    private static long[] getTemperatureLoopTimeDistribution(
        long temperatureLoopTime) throws MPIException {
        LongBuffer mpiOnlyTimingBuffer = ParallelOps.mpiOnlyBuffer;
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.put(temperatureLoopTime);
        ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
        long [] mpiOnlyTimingArray = new long[ParallelOps.worldProcsCount];
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
        return mpiOnlyTimingArray;
    }

    private static void initializeTimers() {
        StressTimings.init(ParallelOps.threadCount);
        StressInternalTimings.init(ParallelOps.threadCount);
        BCInternalTimings.init(ParallelOps.threadCount);
        BCTimings.init(ParallelOps.threadCount);
        MMTimings.init(ParallelOps.threadCount);
    }

    private static void readInitMapping(
        String initialPointsFile, double[][] preX) {
        try (BufferedReader br = Files.newBufferedReader(Paths.get(initialPointsFile),
                                                         Charset.defaultCharset())){
            String line;
            Pattern pattern = Pattern.compile("[\t]");
            int row = 0;
            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line))
                    continue; // continue on empty lines - "while" will break on null anyway;

                String[] splits = pattern.split(line.trim());

                for (int i = 0; i < splits.length; ++i){
                    preX[row][i] = Double.parseDouble(splits[i].trim());
                }
                ++row;
            }
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static String formatElapsedMillis(long elapsed){
        String format = "%dd:%02dH:%02dM:%02dS:%03dmS";
        short millis = (short)(elapsed % (1000.0));
        elapsed = (elapsed - millis) / 1000; // remaining elapsed in seconds
        byte seconds = (byte)(elapsed % 60.0);
        elapsed = (elapsed - seconds) / 60; // remaining elapsed in minutes
        byte minutes =  (byte)(elapsed % 60.0);
        elapsed = (elapsed - minutes) / 60; // remaining elapsed in hours
        byte hours = (byte)(elapsed % 24.0);
        long days = (elapsed - hours) / 24; // remaining elapsed in days
        return String.format(format, days, hours, minutes,  seconds, millis);
    }

    private static void writeOuput(double[][] x, String outputFile)
        throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
        int N = x.length;
        int vecLen = x[0].length;

        DecimalFormat format = new DecimalFormat("#.##########");
        for (int i = 0; i < N; i++) {
            writer.print(String.valueOf(i) + '\t'); // print ID.
            for (int j = 0; j < vecLen; j++) {
                writer.print(format.format(x[i][j]) + '\t'); // print
                // configuration
                // of each axis.
            }
            writer.println("1"); // print label value, which is ONE for all
            // data.
        }
        writer.flush();
        writer.close();

    }

    private static void writeOuput(double[][] X, String labelFile,
                                   String outputFile) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(labelFile));
        String line;
        String parts[];
        Map<String, Integer> labels = new HashMap<>();
        while ((line = reader.readLine()) != null) {
            parts = line.split(" ");
            if (parts.length < 2) {
                // Don't need to throw an error because this is the last part of
                // the computation
                System.out.println("ERROR: Invalid label");
            }
            labels.put(parts[0].trim(), Integer.valueOf(parts[1]));
        }
        reader.close();

        File file = new File(outputFile);
        PrintWriter writer = new PrintWriter(new FileWriter(file));

        int N = X.length;
        int vecLen = X[0].length;

        DecimalFormat format = new DecimalFormat("#.##########");
        for (int i = 0; i < N; i++) {
            writer.print(String.valueOf(i) + '\t'); // print ID.
            for (int j = 0; j < vecLen; j++) {
                writer.print(format.format(X[i][j]) + '\t'); // print
                // configuration
                // of each axis.
            }
            /* TODO Fix bug here - it's from Ryan's code*/
            /*writer.println(labels.get(String.valueOf(ids[i]))); // print label*/
            // value, which
            // is
            // ONE for all data.
        }
        writer.flush();
        writer.close();
    }

    private static void generateVArray(
        short[][] distances, WeightsWrap weights, double[][] vArray) {

        zeroOutArray(vArray);
        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> generateVArrayInternal(
                        threadIdx, distances, weights, vArray[threadIdx])));
        }
        else {
            generateVArrayInternal(0, distances, weights, vArray[0]);
        }
    }

    private static void generateVArrayInternal(
        Integer threadIdx, short[][] distances, WeightsWrap weights, double[] v) {
        int threadRowCount = ParallelOps.threadRowCounts[threadIdx];

        int rowOffset = ParallelOps.threadRowStartOffsets[threadIdx] +
                     ParallelOps.procRowStartOffset;
        for (int i = 0; i < threadRowCount; ++i) {
            int globalRow = i + rowOffset;
            int procLocalRow = globalRow - ParallelOps.procRowStartOffset;
            for (int globalCol = 0; globalCol < ParallelOps.globalColCount; ++globalCol) {
                if (globalRow == globalCol) continue;

                double origD = distances[procLocalRow][globalCol] * INV_SHORT_MAX;
                double weight = weights.getWeight(procLocalRow, globalCol);

                if (origD < 0 || weight == 0) {
                    continue;
                }

                v[i] += weight;
            }
            v[i] += 1;
        }
    }

    private static void calculateConjugateGradient(
        double[][] preX, int targetDimension, int numPoints, double[][] BC, int cgIter, double cgThreshold,
        RefObj<Integer> outCgCount, RefObj<Integer> outRealCGIterations,
        WeightsWrap weights, int blockSize, double[][] vArray, double[][] MMr, double[][] MMAp, double[][][] threadPartialMM)

        throws MPIException {


        zeroOutArray(threadPartialMM);
        CGTimings.startTiming(CGTimings.TimingTask.MM);
        calculateMM(preX, targetDimension, numPoints, weights, blockSize,
                    vArray, MMr, threadPartialMM);

        CGTimings.endTiming(CGTimings.TimingTask.MM);
        // This barrier was necessary for correctness when using
        // a single mmap file
        ParallelOps.worldProcsComm.barrier();

        double[] tmpLHSRow, tmpRHSRow;

        for(int i = 0; i < numPoints; ++i) {
            tmpLHSRow = BC[i];
            tmpRHSRow = MMr[i];
            for (int j = 0; j < targetDimension; ++j) {
                tmpLHSRow[j] -= tmpRHSRow[j];
                tmpRHSRow[j] = tmpLHSRow[j];
            }
        }

        int cgCount = 0;
        CGTimings.startTiming(CGTimings.TimingTask.INNER_PROD);
        double rTr = innerProductCalculation(MMr);
        CGTimings.endTiming(CGTimings.TimingTask.INNER_PROD);
        // Adding relative value test for termination as suggested by Dr. Fox.
        double testEnd = rTr * cgThreshold;

        CGTimings.startTiming(CGTimings.TimingTask.CG_LOOP);
        while(cgCount < cgIter){
            cgCount++;
            outRealCGIterations.setValue(outRealCGIterations.getValue() + 1);

            //calculate alpha
            zeroOutArray(threadPartialMM);
            CGLoopTimings.startTiming(CGLoopTimings.TimingTask.MM);
            calculateMM(BC, targetDimension, numPoints, weights, blockSize,
                        vArray, MMAp, threadPartialMM);
            ParallelOps.worldProcsComm.barrier();
            CGLoopTimings.endTiming(CGLoopTimings.TimingTask.MM);

            CGLoopTimings.startTiming(CGLoopTimings.TimingTask.INNER_PROD_PAP);
            double alpha = rTr
                           /innerProductCalculation(BC, MMAp);
            CGLoopTimings.endTiming(CGLoopTimings.TimingTask.INNER_PROD_PAP);

            //update Xi to Xi+1
            for(int i = 0; i < numPoints; ++i) {
                tmpLHSRow = preX[i];
                tmpRHSRow = BC[i];
                for (int j = 0; j < targetDimension; ++j) {
                    tmpLHSRow[j] += alpha * tmpRHSRow[j];
                }
            }

            if (rTr < testEnd) {
                break;
            }

            //update ri to ri+1
            for(int i = 0; i < numPoints; ++i) {
                tmpLHSRow = MMr[i];
                tmpRHSRow = MMAp[i];
                for (int j = 0; j < targetDimension; ++j) {
                    tmpLHSRow[j] -= alpha * tmpRHSRow[j];
                }
            }

            //calculate beta
            CGLoopTimings.startTiming(CGLoopTimings.TimingTask.INNER_PROD_R);
            double rTr1 = innerProductCalculation(MMr);
            CGLoopTimings.endTiming(CGLoopTimings.TimingTask.INNER_PROD_R);
            double beta = rTr1/rTr;
            rTr = rTr1;

            //update pi to pi+1
            for(int i = 0; i < numPoints; ++i) {
                tmpLHSRow = BC[i];
                tmpRHSRow = MMr[i];
                for (int j = 0; j < targetDimension; ++j) {
                    tmpLHSRow[j] = tmpRHSRow[j] + beta * tmpLHSRow[j];
                }
            }

        }
        CGTimings.endTiming(CGTimings.TimingTask.CG_LOOP);
        outCgCount.setValue(outCgCount.getValue() + cgCount);
    }

    private static void calculateMM(
        double[][] x, int targetDimension, int numPoints, WeightsWrap weights,
        int blockSize, double[][] vArray, double[][] outMM,
        double[][][] internalPartialMM) throws MPIException {

        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> {
                        MMTimings.startTiming(MMTimings.TimingTask.MM_INTERNAL, threadIdx);
                        calculateMMInternal(threadIdx, x, targetDimension,
                                            numPoints, weights, blockSize,
                                            vArray,
                                            internalPartialMM[threadIdx]);
                        MMTimings.endTiming(
                            MMTimings.TimingTask.MM_INTERNAL, threadIdx);
                    }));
        }
        else {
            MMTimings.startTiming(MMTimings.TimingTask.MM_INTERNAL, 0);
            calculateMMInternal(0, x, targetDimension, numPoints, weights,
                                blockSize, vArray, internalPartialMM[0]);
            MMTimings.endTiming(MMTimings.TimingTask.MM_INTERNAL, 0);
        }

        if (ParallelOps.worldProcsCount > 1) {
            MMTimings.startTiming(MMTimings.TimingTask.MM_MERGE, 0);
            mergePartials(internalPartialMM, targetDimension,
                          ParallelOps.mmapXWriteBytes);
            MMTimings.endTiming(MMTimings.TimingTask.MM_MERGE, 0);

            // Important barrier here - as we need to make sure writes are done to the mmap file
            // it's sufficient to wait on ParallelOps.mmapProcComm, but it's cleaner for timings
            // if we wait on the whole world
            ParallelOps.worldProcsComm.barrier();

            if (ParallelOps.isMmapLead) {
                MMTimings.startTiming(MMTimings.TimingTask.COMM, 0);
                ParallelOps.partialXAllGather();
                MMTimings.endTiming(MMTimings.TimingTask.COMM, 0);
            }
            // Each process in a memory group waits here.
            // It's not necessary to wait for a process
            // in another memory map group, hence the use of mmapProcComm.
            // However it's cleaner for any timings to have everyone sync here,
            // so will use worldProcsComm instead.
            ParallelOps.worldProcsComm.barrier();
            MMTimings.startTiming(MMTimings.TimingTask.MM_EXTRACT, 0);
            extractPoints(ParallelOps.fullXBytes,
                                                    ParallelOps.globalColCount,
                                                    targetDimension, outMM);
            MMTimings.endTiming(MMTimings.TimingTask.MM_EXTRACT, 0);
        } else {
            mergePartials(internalPartialMM, targetDimension, outMM);
        }
    }

    private static void calculateMMInternal(
        Integer threadIdx, double[][] x, int targetDimension, int numPoints,
        WeightsWrap weights, int blockSize, double[][] vArray, double[][] outMM) {

        MatrixUtils
            .matrixMultiplyWithThreadOffset(weights, vArray[threadIdx], x,
                ParallelOps.threadRowCounts[threadIdx], targetDimension,
                numPoints, blockSize,
                ParallelOps.threadRowStartOffsets[threadIdx],
                ParallelOps.threadRowStartOffsets[threadIdx]
                + ParallelOps.procRowStartOffset, outMM);
    }


    private static double innerProductCalculation(double[][] a, double[][] b) {
        double sum = 0;
        if (a.length > 0) {
            int row = a.length, col = a[0].length;
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < col; ++j) {
                    sum += a[i][j] * b[i][j];
                }
            }
        }
        return sum;
    }

    private static double innerProductCalculation(double[][] a) {
        double sum = 0;
        if (a.length > 0) {
            int col = a[0].length;
            for (double[] anA : a) {
                for (int j = 0; j < col; ++j) { sum += anA[j] * anA[j]; }
            }
        }
        return sum;
    }

    private static void calculateBC(
        double[][] preX, int targetDimension, double tCur, short[][] distances,
        WeightsWrap weights, int blockSize, double[][] BC,
        double[][][] threadPartialBCInternalBofZ,
        double[][][] threadPartialBCInternalMM)
        throws MPIException, InterruptedException {

        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> {
                        BCTimings.startTiming(BCTimings.TimingTask.BC_INTERNAL,threadIdx);
                        calculateBCInternal(
                            threadIdx, preX, targetDimension, tCur, distances, weights, blockSize, threadPartialBCInternalBofZ[threadIdx], threadPartialBCInternalMM[threadIdx]);
                        BCTimings.endTiming(
                            BCTimings.TimingTask.BC_INTERNAL, threadIdx);
                    }));
        }
        else {
            BCTimings.startTiming(BCTimings.TimingTask.BC_INTERNAL,0);
            calculateBCInternal(
                0, preX, targetDimension, tCur, distances, weights, blockSize,threadPartialBCInternalBofZ[0], threadPartialBCInternalMM[0]);
            BCTimings.endTiming(
                BCTimings.TimingTask.BC_INTERNAL, 0);
        }

        if (ParallelOps.worldProcsCount > 1) {
            BCTimings.startTiming(BCTimings.TimingTask.BC_MERGE, 0);
            mergePartials(threadPartialBCInternalMM, targetDimension,
                          ParallelOps.mmapXWriteBytes);
            BCTimings.endTiming(BCTimings.TimingTask.BC_MERGE, 0);

            // Important barrier here - as we need to make sure writes are done to the mmap file
            // it's sufficient to wait on ParallelOps.mmapProcComm, but it's cleaner for timings
            // if we wait on the whole world
            ParallelOps.worldProcsComm.barrier();

            if (ParallelOps.isMmapLead) {
                BCTimings.startTiming(BCTimings.TimingTask.COMM, 0);
                ParallelOps.partialXAllGather();
                BCTimings.endTiming(BCTimings.TimingTask.COMM, 0);
            }
            // Each process in a memory group waits here.
            // It's not necessary to wait for a process
            // in another memory map group, hence the use of mmapProcComm.
            // However it's cleaner for any timings to have everyone sync here,
            // so will use worldProcsComm instead.
            ParallelOps.worldProcsComm.barrier();

            BCTimings.startTiming(BCTimings.TimingTask.BC_EXTRACT, 0);
            extractPoints(ParallelOps.fullXBytes,
                                                     ParallelOps.globalColCount,
                                                     targetDimension, BC);
            BCTimings.endTiming(BCTimings.TimingTask.BC_EXTRACT, 0);
        } else {
            mergePartials(threadPartialBCInternalMM, targetDimension, BC);
        }
    }

    private static void calculateBCInternal(
        Integer threadIdx, double[][] preX, int targetDimension, double tCur,
        short[][] distances, WeightsWrap weights, int blockSize,
        double[][] internalBofZ, double[][] outMM) {

        BCInternalTimings.startTiming(BCInternalTimings.TimingTask.BOFZ, threadIdx);
        calculateBofZ(threadIdx, preX, targetDimension, tCur,
                                        distances, weights, internalBofZ);
        BCInternalTimings.endTiming(BCInternalTimings.TimingTask.BOFZ, threadIdx);

        // Next we can calculate the BofZ * preX.
        BCInternalTimings.startTiming(BCInternalTimings.TimingTask.MM, threadIdx);
        MatrixUtils.matrixMultiply(internalBofZ, preX, ParallelOps.threadRowCounts[threadIdx],
                                                  targetDimension, ParallelOps.globalColCount, blockSize, outMM);
        BCInternalTimings.endTiming(BCInternalTimings.TimingTask.MM, threadIdx);
    }

    private static void calculateBofZ(
        int threadIdx, double[][] preX, int targetDimension, double tCur, short[][] distances, WeightsWrap weights,
        double[][] outBofZ) {

        int threadRowCount = ParallelOps.threadRowCounts[threadIdx];

        double vBlockValue = -1;

        double diff = 0.0;
        if (tCur > 10E-10) {
            diff = Math.sqrt(2.0 * targetDimension)  * tCur;
        }

        short[] distancesProcLocalRow;
        double[] outBofZLocalRow;
        double origD, weight, dist;
        final int globalRowOffset = ParallelOps.threadRowStartOffsets[threadIdx]
                                    + ParallelOps.procRowStartOffset;
        int globalRow, procLocalRow;
        for (int localRow = 0; localRow < threadRowCount; ++localRow) {
            globalRow = localRow + globalRowOffset;
            procLocalRow = globalRow - ParallelOps.procRowStartOffset;
            outBofZLocalRow = outBofZ[localRow];
            outBofZLocalRow[globalRow] = 0;
            distancesProcLocalRow = distances[procLocalRow];
            for (int globalCol = 0; globalCol < ParallelOps.globalColCount; globalCol++) {
				/*
				 * B_ij = - w_ij * delta_ij / d_ij(Z), if (d_ij(Z) != 0) 0,
				 * otherwise v_ij = - w_ij.
				 *
				 * Therefore, B_ij = v_ij * delta_ij / d_ij(Z). 0 (if d_ij(Z) >=
				 * small threshold) --> the actual meaning is (if d_ij(Z) == 0)
				 * BofZ[i][j] = V[i][j] * deltaMat[i][j] / CalculateDistance(ref
				 * preX, i, j);
				 */
                // this is for the i!=j case. For i==j case will be calculated
                // separately (see above).
                if (globalRow == globalCol) continue;


                origD = distancesProcLocalRow[globalCol] * INV_SHORT_MAX;
                weight = weights.getWeight(procLocalRow,globalCol);

                if (origD < 0 || weight == 0) {
                    continue;
                }

                dist = calculateEuclideanDist(
                    preX[globalRow], preX[globalCol], targetDimension);
                if (dist >= 1.0E-10 && diff < origD) {
                    outBofZLocalRow[globalCol] = (weight * vBlockValue * (origD - diff) / dist);
                } else {
                    outBofZLocalRow[globalCol] = 0;
                }

                outBofZLocalRow[globalRow] -= outBofZLocalRow[globalCol];
            }
        }
    }

    private static void extractPoints(
        Bytes bytes, int numPoints, int dimension, double[][] to) {
        int pos = 0;
        for (int i = 0; i < numPoints; ++i){
            double[] pointsRow = to[i];
            for (int j = 0; j < dimension; ++j) {
                bytes.position(pos);
                pointsRow[j] = bytes.readDouble(pos);
                pos += Double.BYTES;
            }
        }
    }

    private static void mergePartials(double [][][] partials, int dimension, double [][] result){
        int row = 0;
        for (double [][] partial : partials){
            for (double [] point : partial){
                System.arraycopy(point, 0, result[row], 0, dimension);
                ++row;
            }
        }
    }

    private static void mergePartials(
        double[][][] partials, int targetDimension, Bytes result){
        result.position(0);
        for (double [][] partial : partials){
            for (double [] point : partial){
                for (int i = 0; i < targetDimension; ++i){
                    result.writeDouble(point[i]);
                }
            }
        }
    }

    private static double calculateStress(
        double[][] preX, int targetDimension, double tCur, short[][] distances,
        WeightsWrap weights, double invSumOfSquareDist, double[] internalPartialSigma)
        throws MPIException {

        double stress = 0.0;

        if (ParallelOps.threadCount > 1) {
            IntStream.range(0, ParallelOps.threadCount).forEach(i -> internalPartialSigma[i] = 0.0);
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> {
                        StressTimings.startTiming(StressTimings.TimingTask.STRESS_INTERNAL, threadIdx);
                        internalPartialSigma[threadIdx] =
                        calculateStressInternal(threadIdx, preX, targetDimension, tCur,

                                                distances, weights);
                        StressTimings.endTiming(StressTimings.TimingTask.STRESS_INTERNAL, threadIdx);
                    }));
            // Sum across threads and accumulate to stress
            for (int i = 0; i < ParallelOps.threadCount; ++i){
                stress += internalPartialSigma[i];
            }
        }
        else {
            stress = calculateStressInternal(0, preX, targetDimension, tCur,
                                                     distances, weights);
        }

        if (ParallelOps.worldProcsCount > 1) {
            /*
             *//*Write thread local reduction to shared memory map*//*
            ParallelOps.mmapSWriteBytes.position(0);
            ParallelOps.mmapSWriteBytes.writeDouble(stress);

            // Important barrier here - as we need to make sure writes are done
            // to the mmap file.
            // It's sufficient to wait on ParallelOps.mmapProcComm,
            // but it's cleaner for timings if we wait on the whole world
            ParallelOps.worldProcsComm.barrier();
            if (ParallelOps.isMmapLead) {
                 *//*Node local reduction using shared memory maps*//*
                ParallelOps.mmapSReadBytes.position(0);
                stress = 0.0;
                for (int i = 0; i < ParallelOps.mmapProcsCount; ++i){
                    stress += ParallelOps.mmapSReadBytes.readDouble();
                }
                ParallelOps.mmapSWriteBytes.position(0);
                ParallelOps.mmapSWriteBytes.writeDouble(stress);

                 *//*Leaders participate in MPI AllReduce*//*
                StressTimings.startTiming(StressTimings.TimingTask.COMM, 0);
                ParallelOps.partialSAllReduce(MPI.SUM);
                StressTimings.endTiming(StressTimings.TimingTask.COMM, 0);
            }

            // Each process in a memory group waits here.
            // It's not necessary to wait for a process
            // in another memory map group, hence the use of mmapProcComm.
            // However it's cleaner for any timings to have everyone sync here,
            // so will use worldProcsComm instead.
            ParallelOps.worldProcsComm.barrier();
            ParallelOps.mmapSReadBytes.position(0);
            stress = ParallelOps.mmapSReadBytes.readDouble();
*/
            /* Note - This works!, so have to see why mmap version failed*/
            // Barrier for cleaner timings
            ParallelOps.worldProcsComm.barrier();
            StressTimings.startTiming(StressTimings.TimingTask.COMM, 0);
            stress = ParallelOps.allReduce(stress);
            StressTimings.endTiming(StressTimings.TimingTask.COMM, 0);
        }
        return stress * invSumOfSquareDist;
    }

    private static double calculateStressInternal(
        int threadIdx, double[][] preX, int targetDim, double tCur, short[][] distances, WeightsWrap weights) {

        StressInternalTimings.startTiming(StressInternalTimings.TimingTask.COMP, threadIdx);
        double sigma = 0.0;
        double diff = 0.0;
        if (tCur > 10E-10) {
            diff = Math.sqrt(2.0 * targetDim) * tCur;
        }

        int threadRowCount = ParallelOps.threadRowCounts[threadIdx];
        final int globalRowOffset = ParallelOps.threadRowStartOffsets[threadIdx]
                                    + ParallelOps.procRowStartOffset;

        int globalRow, procLocalRow;
        short[] distancesProcLocalRow;
        double origD, weight, euclideanD;
        double heatD, tmpD;
        for (int localRow = 0; localRow < threadRowCount; ++localRow){
            globalRow = localRow + globalRowOffset;
            procLocalRow = globalRow - ParallelOps.procRowStartOffset;
            distancesProcLocalRow = distances[procLocalRow];
            for (int globalCol = 0; globalCol < ParallelOps.globalColCount; globalCol++) {
                origD = distancesProcLocalRow[globalCol] * INV_SHORT_MAX;
                weight = weights.getWeight(procLocalRow,globalCol);

                if (origD < 0 || weight == 0) {
                    continue;
                }

                euclideanD = globalRow != globalCol ? calculateEuclideanDist(
                    preX[globalRow], preX[globalCol], targetDim) : 0.0;

                heatD = origD - diff;
                tmpD = origD >= diff ? heatD - euclideanD : -euclideanD;
                sigma += weight * tmpD * tmpD;
            }
        }
        StressInternalTimings.endTiming(StressInternalTimings.TimingTask.COMP, threadIdx);
        return sigma;
    }

    private static double calculateEuclideanDist(
        double[] v, double[] w, int targetDim) {
        double dist = 0.0;
        for (int k = 0; k < targetDim; k++) {
            try {
                double diff = v[k] - w[k];
                dist += diff * diff;
            } catch (IndexOutOfBoundsException e){
                // Usually this should not happen, also this is not
                // necessary to catch, but it was found that some parent block
                // in HJ hides this error if/when it happens, so explicitly
                // printing it here.
                e.printStackTrace();
            }
        }

        dist = Math.sqrt(dist);
        return dist;
    }

    private static void generateInitMapping(
        int numPoints, int targetDim, double[][] preX) throws MPIException {

        Bytes fullBytes = ParallelOps.fullXBytes;
        if (ParallelOps.worldProcRank == 0) {
            int pos = 0;
            // Use Random class for generating random initial mapping solution.
            Random rand = new Random(System.currentTimeMillis());
            for (int i = 0; i < numPoints; i++) {
                for (int j = 0; j < targetDim; j++) {
                    fullBytes.position(pos);
                    fullBytes.writeDouble(rand.nextBoolean()
                                              ? rand.nextDouble()
                                              : -rand.nextDouble());
                    pos += Double.BYTES;
                }
            }
        }

        if (ParallelOps.worldProcsCount > 1){
            // Broadcast initial mapping to others
            ParallelOps.broadcast(ParallelOps.fullXByteBuffer, numPoints * targetDim*Double.BYTES, 0);
        }
        extractPoints(fullBytes, numPoints, targetDim, preX);
    }

    private static DoubleStatistics calculateStatistics(
        short[][] distances, WeightsWrap weights, RefObj<Integer> missingDistCount) throws MPIException {

        final DoubleStatistics[] threadDistanceSummaries =
            new DoubleStatistics[ParallelOps.threadCount];
        final int [] missingDistCounts = new int[ParallelOps.threadCount];
        IntStream.range(0, ParallelOps.threadCount).forEach(i -> missingDistCounts[i] = 0);

        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> threadDistanceSummaries[threadIdx] =
                        calculateStatisticsInternal(
                            threadIdx, distances, weights, missingDistCounts)));

            // Sum across threads and accumulate to zeroth entry
            IntStream.range(1, ParallelOps.threadCount).forEach(
                i -> {
                    threadDistanceSummaries[0]
                        .combine(threadDistanceSummaries[i]);
                    missingDistCounts[0] += missingDistCounts[i];
                });
        }
        else {
            threadDistanceSummaries[0] = calculateStatisticsInternal(
                0, distances, weights, missingDistCounts);
        }

        if (ParallelOps.worldProcsCount > 1) {
            threadDistanceSummaries[0] =
                ParallelOps.allReduce(threadDistanceSummaries[0]);
            missingDistCounts[0] = ParallelOps.allReduce(missingDistCounts[0]);
        }
        missingDistCount.setValue(missingDistCounts[0]);
        return threadDistanceSummaries[0];
    }

    private static TransformationFunction loadFunction(String classFile) {
        ClassLoader classLoader = Program.class.getClassLoader();
        try {
            Class aClass = classLoader.loadClass(classFile);
            return (TransformationFunction) aClass.newInstance();
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException e) {
            throw new RuntimeException("Failed to load class: " + classFile, e);
        }
    }

    private static void readDistancesAndWeights(boolean isSammon) {
        TransformationFunction function;
        if (!Strings.isNullOrEmpty(config.transformationFunction)) {
            function = loadFunction(config.transformationFunction);
        } else {
            function = (config.distanceTransform != 1.0
                            ? (d -> Math.pow(d, config.distanceTransform))
                            : null);
        }

        distances = config.repetitions == 1 ? BinaryReader2D.readRowRange(
            config.distanceMatrixFile, ParallelOps.procRowRange,
            ParallelOps.globalColCount, byteOrder, true, function)
            : BinaryReader2D.readRowRange(
                config.distanceMatrixFile, ParallelOps.procRowRange,
                ParallelOps.globalColCount, byteOrder, true, function, config.repetitions);

        short[][] w = null;
        if (!Strings.isNullOrEmpty(config.weightMatrixFile)){
            function = !Strings.isNullOrEmpty(config.weightTransformationFunction)
                ? loadFunction(config.weightTransformationFunction)
                : null;
            w = config.repetitions == 1 ? BinaryReader2D
                .readRowRange(config.weightMatrixFile, ParallelOps.procRowRange,
                              ParallelOps.globalColCount, byteOrder, true,
                              function)
            : BinaryReader2D
                .readRowRange(config.weightMatrixFile, ParallelOps.procRowRange,
                    ParallelOps.globalColCount, byteOrder, true,
                    function, config.repetitions);
        }
        weights = new WeightsWrap(w, distances, isSammon);
    }

    private static DoubleStatistics calculateStatisticsInternal(
        int threadIdx, short[][] distances, WeightsWrap weights, int[] missingDistCounts) {

        DoubleStatistics stat = new DoubleStatistics();

        int threadRowCount = ParallelOps.threadRowCounts[threadIdx];
        final int threadRowStartOffset = ParallelOps.threadRowStartOffsets[threadIdx];

        int procLocalRow;
        short[] distancesProcLocalRow;
        double origD, weight;
        for (int localRow = 0; localRow < threadRowCount; ++localRow){
            procLocalRow = localRow + threadRowStartOffset;
            distancesProcLocalRow = distances[procLocalRow];
            for (int globalCol = 0; globalCol < ParallelOps.globalColCount; globalCol++) {
                origD = distancesProcLocalRow[globalCol] * INV_SHORT_MAX;
                weight = weights.getWeight(procLocalRow,globalCol);
                if (origD < 0) {
                    // Missing distance
                    ++missingDistCounts[threadIdx];
                    continue;
                }
                if (weight == 0) continue; // Ignore zero weights

                stat.accept(origD);
            }
        }

        return stat;
    }

    private static void readConfiguration(CommandLine cmd) {
        config = ConfigurationMgr.LoadConfiguration(
            cmd.getOptionValue(Constants.CMD_OPTION_LONG_C)).damdsSection;
        ParallelOps.nodeCount =
            Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_N));
        ParallelOps.threadCount =
            Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_T));
        ParallelOps.mmapsPerNode = cmd.hasOption(Constants.CMD_OPTION_SHORT_MMAPS) ? Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_SHORT_MMAPS)) : 1;
        ParallelOps.mmapScratchDir = cmd.hasOption(Constants.CMD_OPTION_SHORT_MMAP_SCRATCH_DIR) ? cmd.getOptionValue(Constants.CMD_OPTION_SHORT_MMAP_SCRATCH_DIR) : ".";

        byteOrder =
            config.isBigEndian ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN;
        BlockSize = config.blockSize;
    }

    /**
     * Parse command line arguments
     *
     * @param args Command line arguments
     * @param opts Command line options
     * @return An <code>Optional&lt;CommandLine&gt;</code> object
     */
    private static Optional<CommandLine> parseCommandLineArguments(
        String[] args, Options opts) {

        CommandLineParser optParser = new GnuParser();

        try {
            return Optional.fromNullable(optParser.parse(opts, args));
        }
        catch (ParseException e) {
            e.printStackTrace();
        }
        return Optional.fromNullable(null);
    }
}
