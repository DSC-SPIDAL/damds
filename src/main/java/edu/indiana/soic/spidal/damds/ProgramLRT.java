package edu.indiana.soic.spidal.damds;

import com.google.common.base.Optional;
import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.configuration.ConfigurationMgr;
import edu.indiana.soic.spidal.configuration.section.DAMDSSection;
import edu.indiana.soic.spidal.damds.threads.SpidalThreads;
import edu.indiana.soic.spidal.damds.threads.ThreadBitAssigner;
import mpi.MPIException;
import net.openhft.affinity.Affinity;
import org.apache.commons.cli.*;

import java.io.*;
import java.nio.ByteOrder;
import java.nio.LongBuffer;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;


public class ProgramLRT {
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
        programOptions.addOption(
                Constants.CMD_OPTION_SHORT_BIND_THREADS, true,
                Constants.CMD_OPTION_DESCRIPTION_BIND_THREADS);

        programOptions.addOption(
                Constants.CMD_OPTION_SHORT_CPS, true,
                Constants.CMD_OPTION_DESCRIPTION_CPS);


    }

    //Config Settings
    public static DAMDSSection config;
    public static ByteOrder byteOrder;

    public static int BlockSize;
    private static Utils utils = new Utils(0);
    private static boolean bind;
    private static int cps;

    private static SpidalThreads threads = null;

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
            ParallelOps.setParallelDecomposition(
                config.numberDataPoints, config.targetDimension);

            /*if (ParallelOps.threadCount > 1) {
                threads = new SpidalThreads(ParallelOps.threadCount, false, true,
                        48, ParallelOps.worldProcRank * 12 + 1);
            }*/
            // Note - a barrier to get cleaner timings
            ParallelOps.worldProcsComm.barrier();
            Stopwatch mainTimer = Stopwatch.createStarted();

            utils.printMessage("\n== DAMDS run started on " + new Date() + " ==\n");
            utils.printMessage(config.toString(false));

            /* TODO - Fork - join starts here */
            if (ParallelOps.threadCount > 1) {
                Lock lock = new ReentrantLock();
                launchHabaneroApp(
                    () -> forallChunked(
                        0, ParallelOps.threadCount - 1,
                        (threadIdx) -> {

                            if (bind) {
                                BitSet bitSet = ThreadBitAssigner.getBitSet(ParallelOps.worldProcRank, threadIdx, ParallelOps.threadCount, cps);
                                Affinity.setAffinity(bitSet);
                            }

                            final ProgramWorker worker = new ProgramWorker
                                    (threadIdx,
                                    ParallelOps
                                    .threadComm, config, byteOrder,
                                            BlockSize, mainTimer,
                                            lock);
                            try {
                                worker.run();
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }));
                /*threads.forall(
                        (threadIdx) -> {
                            new ProgramWorker(threadIdx, ParallelOps
                                    .threadComm, config, byteOrder,
                                    BlockSize, mainTimer,lock).run();
                        });*/
            }
            else {
                if (bind) {
                    BitSet bitSet = ThreadBitAssigner.getBitSet(ParallelOps.worldProcRank, 0, ParallelOps.threadCount, cps);
                    Affinity.setAffinity(bitSet);
                }
                new ProgramWorker(0, ParallelOps.threadComm, config,
                        byteOrder, BlockSize, mainTimer, null)
                        .run();
            }


            /* TODO - Fork-join should end here */


            /*
             // TODO Fix error handling here
            if (Strings.isNullOrEmpty(config.labelFile) || config.labelFile.toUpperCase().endsWith(
                "NOLABEL")) {
                try {
                    writeOuput(preX, config.targetDimension, config.pointsFile);
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
            } else {
                try {
                    writeOuput(preX, config.targetDimension, config.labelFile, config.pointsFile);
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
            }

            Double finalStress = calculateStress(
                preX, config.targetDimension, tCur, distances, weights,
                INV_SUM_OF_SQUARE, partialSigma);

            mainTimer.stop();


            utils.printMessage("Finishing DAMDS run ...");
            long totalTime = mainTimer.elapsed(TimeUnit.MILLISECONDS);
            long temperatureLoopTime = loopTimer.elapsed(TimeUnit.MILLISECONDS);
            utils.printMessage(
                String.format(
                    "  Total Time: %s (%d ms) Loop Time: %s (%d ms)",
                    formatElapsedMillis(totalTime), totalTime,
                    formatElapsedMillis(temperatureLoopTime), temperatureLoopTime));
            utils.printMessage("  Total Loops: " + loopNum);
            utils.printMessage("  Total Iterations: " + smacofRealIterations);
            utils.printMessage(
                String.format(
                    "  Total CG Iterations: %d Avg. CG Iterations: %.5g",
                    outRealCGIterations.getValue(), (outRealCGIterations.getValue() * 1.0) / smacofRealIterations));
            utils.printMessage("  Final Stress:\t" + finalStress);

            printTimings(totalTime, temperatureLoopTime);*/

            utils.printMessage("== DAMDS run completed on " + new Date() + " ==");

            /*if (ParallelOps.threadCount > 1) {
                threads.shutDown();
            }*/
            ParallelOps.tearDownParallelism();
        }
        catch (MPIException | IOException e) {
            utils.printAndThrowRuntimeException(new RuntimeException(e));
        }
    }

    /*private static void printTimings(long totalTime, long temperatureLoopTime)
        throws MPIException {
        String mainHeader =
            "  Total\tTempLoop\tPreStress\tStressLoop\tBC\tCG\tStress" +
            "\tBCInternal\tBComm\tBCInternalBofZ\tBCInternalMM\tCGMM" +
            "\tCGInnerProd\tCGLoop\tCGLoopMM\tCGLoopInnerProdPAP" +
            "\tCGLoopInnerProdR\tMMInternal\tMMComm\tBCMerge\tBCExtract\tMMMerge\tMMExtract\tStressInternal\tStressComm\tStressInternalComp";
        utils.printMessage(
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

        utils.printMessage(
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
        utils.printMessage(
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
        utils.printMessage(
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
    }*/

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

    private static void writeOuput(double[] x, int vecLen, String outputFile)
        throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
        int N = x.length / vecLen;

        DecimalFormat format = new DecimalFormat("#.##########");
        for (int i = 0; i < N; i++) {
            int index = i * vecLen;
            writer.print(String.valueOf(i) + '\t'); // print ID.
            for (int j = 0; j < vecLen; j++) {
                writer.print(format.format(x[index + j]) + '\t'); // print
                // configuration
                // of each axis.
            }
            writer.println("1"); // print label value, which is ONE for all
            // data.
        }
        writer.flush();
        writer.close();

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

        bind = !cmd.hasOption(Constants.CMD_OPTION_SHORT_BIND_THREADS) ||
                Boolean.parseBoolean(cmd.getOptionValue(Constants.CMD_OPTION_SHORT_BIND_THREADS));
        cps = (cmd.hasOption(Constants.CMD_OPTION_SHORT_CPS)) ? Integer.parseInt(Constants.CMD_OPTION_SHORT_CPS) : -1;
        if (cps == -1){
            utils.printMessage("Disabling thread binding as cps is not specified");
            bind = false;
        }
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
