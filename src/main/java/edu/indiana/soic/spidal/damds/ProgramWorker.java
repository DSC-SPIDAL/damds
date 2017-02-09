package edu.indiana.soic.spidal.damds;

import com.google.common.base.Optional;
import com.google.common.base.Stopwatch;
import com.google.common.base.Strings;
import edu.indiana.soic.spidal.common.*;
import edu.indiana.soic.spidal.configuration.section.DAMDSSection;
import edu.indiana.soic.spidal.damds.threads.ThreadCommunicator;
import net.openhft.lang.io.Bytes;
import org.apache.commons.cli.*;

import java.io.*;
import java.nio.ByteOrder;
import java.nio.LongBuffer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.regex.Pattern;

public class ProgramWorker {
    // Constants
    private final double INV_SHORT_MAX = 1.0 / Short.MAX_VALUE;
    private final double SHORT_MAX = Short.MAX_VALUE;

    // Calculated Constants
    private double INV_SUM_OF_SQUARE;

    // Arrays
    private double[] preX;
    private double[] BC;
    private double[] MMr;
    private double[] MMAp;

    private double[][] threadPartialBofZ;
    private double[] threadPartialMM;

    private double[] v;

    //Config Settings
    private DAMDSSection config;
    private ByteOrder byteOrder;
    private short[] distances;
    private WeightsWrap1D weights;

    private int BlockSize;

    private int threadId;
    private Range globalThreadRowRange;
    private Range threadLocalRowRange;

    final private RefObj<Integer> refInt = new RefObj<>();
    final private RefObj<Double> refDouble = new RefObj<>();

    private ThreadCommunicator threadComm;
    private Utils utils;
    private Stopwatch mainTimer;


    private Bytes threadLocalFullXBytes;
    private Bytes threadLocalMmapXWriteBytes;

    private Lock lock;

    public ProgramWorker(int threadId, ThreadCommunicator comm, DAMDSSection
            config, ByteOrder byteOrder, int blockSize, Stopwatch mainTimer,
                         Lock lock) {
        this.threadId = threadId;
        this.threadComm = comm;
        this.config = config;
        this.byteOrder = byteOrder;
        this.BlockSize = blockSize;
        this.mainTimer = mainTimer;
        this.lock = lock;
        utils = new Utils(threadId);

    }

    public void setup() {
        final int threadRowCount = config.numberDataPoints;
        final int threadLocalRowStartOffset =0;
        final int globalThreadRowStartOffset = 0;
        globalThreadRowRange = new Range(
                globalThreadRowStartOffset, globalThreadRowStartOffset +
                threadRowCount - 1);
        threadLocalRowRange = new Range(
                threadLocalRowStartOffset, (
                threadLocalRowStartOffset + threadRowCount - 1));

        if (lock != null) {
            lock.lock();
        }




        if (lock != null) {
            lock.unlock();
        }
    }

    public void run() throws IOException {
        try {
            setup();
            readDistancesAndWeights(config.isSammon);
            /*System.out.println("Rank " + ParallelOps.worldProcRank + " " +
                    "TID " + threadId + "Came " +
                    "here ");*/

            RefObj<Integer> missingDistCount = new RefObj<>();


            // Allocating point arrays once for all
            allocateArrays();

            if (Strings.isNullOrEmpty(config.initialPointsFile)) {
                System.out.println("***********ERROR initial point file null");
            } else {
                readInitMapping(config.initialPointsFile, preX, config
                        .targetDimension);
            }


            double tCur = 0.0;
            double tMax = 1.0 / Math.sqrt(2.0 * config
                    .targetDimension);
            double tMin = config.tMinFactor * 1.0 / Math.sqrt(2.0 * config.targetDimension);

            generateV(distances, weights, v);

            utils.printMessage("\nInitial stress=" );

            tCur = config.alpha * tMax;

            // TODO - turning off barriers
            /*if (threadId == 0) {
                ParallelOps.worldProcsComm.barrier();
            }
            threadComm.barrier();*/
            if (threadId == 0) {
                mainTimer.stop();
                utils.printMessage("\nUp to the loop took " + mainTimer.elapsed(
                        TimeUnit.SECONDS) + " seconds");
                mainTimer.start();
            }

            Stopwatch loopTimer = Stopwatch.createStarted();

            int loopNum = 0;
            double diffStress;
            double stress = -1.0;
            RefObj<Integer> outRealCGIterations = new RefObj<>(0);
            RefObj<Integer> cgCount = new RefObj<>(0);
            int smacofRealIterations = 0;
            while (true) {



                diffStress = config.threshold + 1.0;

                utils.printMessage(
                        String.format(
                                "\nStart of loop %d Temperature (T_Cur) %.5g",
                                loopNum, tCur));

                int itrNum = 0;
                cgCount.setValue(0);
                while (diffStress >= config.threshold) {

                    zeroOutArray(threadPartialMM);
                    calculateBC(
                            preX, config.targetDimension, tCur, distances,
                            weights, BlockSize, BC, threadPartialBofZ,
                            threadPartialMM);

                    // TODO - turning off barriers
                    /*if (threadId == 0) {
                        // This barrier was necessary for correctness when using
                        // a single mmap file
                        ParallelOps.worldProcsComm.barrier();
                    }
                    threadComm.barrier();*/

                    /*stressLoopTimings.startTiming(
                            StressLoopTimings.TimingTask.CG);
                    calculateConjugateGradient(preX, config.targetDimension,
                            config.numberDataPoints,
                            BC,
                            config.cgIter,
                            config.cgErrorThreshold, cgCount,
                            outRealCGIterations, weights,
                            BlockSize, v, MMr, MMAp, threadPartialMM);
                    stressLoopTimings.endTiming(
                            StressLoopTimings.TimingTask.CG);*/


                    /*stressLoopTimings.startTiming(
                            StressLoopTimings.TimingTask.STRESS);
                    stress = calculateStress(
                            preX, config.targetDimension, tCur, distances,
                            weights,
                            INV_SUM_OF_SQUARE);
                    stressLoopTimings.endTiming(
                            StressLoopTimings.TimingTask.STRESS);*/

                    diffStress = 1.0 - stress;
                    diffStress = stress;

                    if ((itrNum % 10 == 0) || (itrNum >= config.stressIter)) {
                        utils.printMessage(
                                String.format(
                                        "  Loop %d Iteration %d Avg CG count " +
                                                "%.5g " +
                                                "Stress " +
                                                "%.5g", loopNum, itrNum,
                                        (cgCount.getValue() * 1.0 / (itrNum +
                                                1)),
                                        stress));
                    }
                    ++itrNum;
                    ++smacofRealIterations;
                }

                --itrNum;
                if (itrNum >= 0 && !(itrNum % 10 == 0) && !(itrNum >=
                        config.stressIter)) {
                    utils.printMessage(
                            String.format(
                                    "  Loop %d Iteration %d Avg CG count %.5g" +
                                            " " +
                                            "Stress %.5g",
                                    loopNum, itrNum,
                                    (cgCount.getValue() * 1.0 / (itrNum + 1))
                                    , stress));
                }

                utils.printMessage(
                        String.format(
                                "End of loop %d Total Iterations %d Avg CG " +
                                        "count %.5g" +
                                        " Stress %.5g",
                                loopNum, (itrNum + 1),
                                (cgCount.getValue() * 1.0 / (itrNum + 1)),
                                stress));

                if (tCur == 0)
                    break;
                tCur *= config.alpha;
                if (tCur < tMin)
                    tCur = 0;
                ++loopNum;

                /* Note - quick way to test programs without running full
                * number of temperature loops */
                if (config.maxtemploops > 0 && loopNum == config.maxtemploops) {
                    break;
                }
            }
            loopTimer.stop();

            double QoR1 = stress / (config.numberDataPoints * (config
                    .numberDataPoints - 1) / 2);
            double QoR2 = 1.0;

            utils.printMessage(
                    String.format(
                            "Normalize1 = %.5g Normalize2 = %.5g", QoR1, QoR2));
            utils.printMessage(
                    String.format(
                            "Average of Delta(original distance) = %.5g"));



            if (threadId == 0) {
                if (0 == 0){
                    Utils.writeOutput(preX, ProgramLRT.config.targetDimension, ProgramLRT.config.pointsFile);
                }
                mainTimer.stop();
            }

            utils.printMessage("Finishing DAMDS run ...");
            long totalTime = mainTimer.elapsed(TimeUnit.MILLISECONDS);
            long temperatureLoopTime = loopTimer.elapsed(TimeUnit.MILLISECONDS);
            utils.printMessage(
                    String.format(
                            "  Total Time: %s (%d ms) Loop Time: %s (%d ms)",
                            formatElapsedMillis(totalTime), totalTime,
                            formatElapsedMillis(temperatureLoopTime),
                            temperatureLoopTime));
            utils.printMessage("  Total Loops: " + loopNum);
            utils.printMessage("  Total Iterations: " + smacofRealIterations);
            utils.printMessage(
                    String.format(
                            "  Total CG Iterations: %d Avg. CG Iterations: %" +
                                    ".5g",
                            outRealCGIterations.getValue(),
                            (outRealCGIterations.getValue() * 1.0) /
                                    smacofRealIterations));
            utils.printMessage("  Final Stress:\t" );
            // TODO - fix print timings
            /*printTimings(totalTime, temperatureLoopTime);*/
            threadComm.barrier();
        } catch (Exception e) {
            utils.printAndThrowRuntimeException(new RuntimeException(e));
        }
    }



    private static void prettyPrintArray(
            String title, double[] vals, PrintWriter printWriter) {
        String str;
        printWriter.println(title);
        str = Arrays.toString(vals);
        printWriter.println(
                str.substring(1, str.length() - 1).replace(',', '\t'));
        printWriter.println();
    }

    private void allocateArrays() {
        // Allocating point arrays once for all
        final int numberDataPoints = config.numberDataPoints;
        final int targetDimension = config.targetDimension;

        preX = new double[numberDataPoints * targetDimension];
        BC = new double[numberDataPoints * targetDimension];
        MMr = new double[numberDataPoints * targetDimension];
        MMAp = new double[numberDataPoints * targetDimension];
        final int threadRowCount = config.numberDataPoints;
        threadPartialBofZ = new double[threadRowCount][config.numberDataPoints];
        threadPartialMM = new double[threadRowCount
                * config.targetDimension];
        v = new double[threadRowCount];
    }

    private void zeroOutArray(double[] a) {
        Arrays.fill(a, 0.0d);
    }



    private void changeZeroDistancesToPostiveMin(
            short[] distances, double positiveMin) {
        double tmpD;
        for (int i = 0; i < distances.length; ++i) {
            tmpD = distances[i] * INV_SHORT_MAX;
            if (tmpD < positiveMin && tmpD >= 0.0) {
                distances[i] = (short) (positiveMin * SHORT_MAX);
            }
        }
    }



    private void readInitMapping(
            String initialPointsFile, double[] preX, int dimension)
            throws BrokenBarrierException, InterruptedException {
        if (threadId == 0) {
            try (BufferedReader br = Files
                    .newBufferedReader(Paths.get(initialPointsFile),
                            Charset.defaultCharset())) {
                String line;
                Pattern pattern = Pattern.compile("[\t]");
                int row = 0;
                while ((line = br.readLine()) != null) {

                    if (Strings.isNullOrEmpty(line)) {
                        continue; // continue on empty lines - "while" will
                    }
                    // break on null anyway;

                    String[] splits = pattern.split(line.trim());

                    for (int i = 0; i < dimension; ++i) {
                        preX[row + i] = Double.parseDouble(splits[i].trim());
                    }
                    row += dimension;

                }
            } catch (Exception e) {
                throw new RuntimeException(e);
            }

        }
        threadComm.barrier();
        threadComm.bcastDoubleArrayOverThreads(threadId, preX, 0);
    }

    public static String formatElapsedMillis(long elapsed) {
        String format = "%dd:%02dH:%02dM:%02dS:%03dmS";
        short millis = (short) (elapsed % (1000.0));
        elapsed = (elapsed - millis) / 1000; // remaining elapsed in seconds
        byte seconds = (byte) (elapsed % 60.0);
        elapsed = (elapsed - seconds) / 60; // remaining elapsed in minutes
        byte minutes = (byte) (elapsed % 60.0);
        elapsed = (elapsed - minutes) / 60; // remaining elapsed in hours
        byte hours = (byte) (elapsed % 24.0);
        long days = (elapsed - hours) / 24; // remaining elapsed in days
        return String.format(format, days, hours, minutes, seconds, millis);
    }



    private static void writeOuput(double[] X, int vecLen, String labelFile,
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

        int N = X.length / 3;

        DecimalFormat format = new DecimalFormat("#.##########");
        for (int i = 0; i < N; i++) {
            int index = i * vecLen;
            writer.print(String.valueOf(i) + '\t'); // print ID.
            for (int j = 0; j < vecLen; j++) {
                writer.print(format.format(X[index + j]) + '\t'); // print
                // configuration
                // of each axis.
            }
            /* TODO Fix bug here - it's from Ryan's code*/
            /*writer.println(labels.get(String.valueOf(ids[i]))); // print
            label*/
            // value, which
            // is
            // ONE for all data.
        }
        writer.flush();
        writer.close();
    }

    private void generateV(
            short[] distances, WeightsWrap1D weights, double[] v) {
        zeroOutArray(v);
        int threadRowCount = config.numberDataPoints;

        int rowOffset = 0;
        for (int threadLocalRow = 0; threadLocalRow < threadRowCount;
             ++threadLocalRow) {
            int globalRow = threadLocalRow + rowOffset;
            for (int globalCol = 0; globalCol < config.numberDataPoints;
                 ++globalCol) {
                if (globalRow == globalCol) continue;

                double origD = distances[threadLocalRow * config.numberDataPoints + globalCol] * INV_SHORT_MAX;
//                double weight = weights.getWeight(threadLocalRow, globalCol);
                double weight = 1.0;

                if (origD < 0 || weight == 0) {
                    continue;
                }

                v[threadLocalRow] += weight;
            }
            v[threadLocalRow] += 1;
        }

    }






    private void calculateMMInternal(
            double[] x, int targetDimension, int numPoints,
            WeightsWrap1D weights, int blockSize, double[] v, double[] outMM) {

        MatrixUtils.matrixMultiplyWithThreadOffset(weights, v, x,
                globalThreadRowRange.getLength(), targetDimension, numPoints,
                blockSize, 0, globalThreadRowRange.getStartIndex(), outMM);
    }


    private static double innerProductCalculation(double[] a, double[] b) {
        double sum = 0;
        if (a.length > 0) {
            for (int i = 0; i < a.length; ++i) {
                sum += a[i] * b[i];
            }
        }
        return sum;
    }

    private static double innerProductCalculation(double[] a) {
        double sum = 0.0;
        if (a.length > 0) {
            for (double anA : a) {
                sum += anA * anA;
            }
        }
        return sum;
    }

    private void calculateBC(
            double[] preX, int targetDimension, double tCur, short[] distances,
            WeightsWrap1D weights, int blockSize, double[] BC,
            double[][] threadPartialBCInternalBofZ,
            double[] threadPartialBCInternalMM)
            throws Exception, InterruptedException, BrokenBarrierException {

        // TODO: VT DEBUG
        calculateBCInternal(
                preX, targetDimension, tCur, distances, weights, blockSize,
                threadPartialBCInternalBofZ, threadPartialBCInternalMM);



        //System.out.println(threadId);
        // TODO: VT DEBUG
        /*threadComm.collect2(0,
                threadPartialBCInternalMM, threadLocalMmapXWriteBytes, threadId);*/

        // TODO: VT DEBUG
//        if (ParallelOps.worldProcsCount > 1) {
//            if (threadId == 0) {
//                /*bcTimings.startTiming(BCTimings.TimingTask.COMM);
//                ParallelOps.allGather();
//                bcTimings.endTiming(BCTimings.TimingTask.COMM);*/
//
//                // Important barrier here - as we need to make sure writes
//                // are done to the mmap file
//
//                // it's sufficient to wait on ParallelOps.mmapProcComm, but
//                // it's cleaner for timings
//                // if we wait on the whole world
//                ParallelOps.worldProcsComm.barrier();
//
//                if (ParallelOps.isMmapLead) {
//                    // bcTimings.startTiming(BCTimings.TimingTask.COMM);
//                    ParallelOps.partialXAllGather();
//                    // bcTimings.endTiming(BCTimings.TimingTask.COMM, 0);
//                }
//                // Each process in a memory group waits here.
//                // It's not necessary to wait for a process
//                // in another memory map group, hence the use of
//                // mmapProcComm.
//                // However it's cleaner for any timings to have everyone sync
//                // here, so will use worldProcsComm instead.
//                ParallelOps.worldProcsComm.barrier();
//            }
//            threadComm.barrier();
//        }
        // TODO: VT DEBUG
        /*threadComm.copy2(ParallelOps.worldProcsCount > 1
                        ? threadLocalFullXBytes
                        : threadLocalMmapXWriteBytes, BC,
                ParallelOps.globalColCount * targetDimension, threadId);*/
        //threadComm.barrier();


    }

    private void calculateBCInternal(
            double[] preX, int targetDimension, double tCur, short[] distances,
            WeightsWrap1D weights, int blockSize, double[][] internalBofZ,
            double[] outMM) {

        calculateBofZ(preX, targetDimension, tCur,
                distances, weights, internalBofZ);

        // Next we can calculate the BofZ * preX.
        // TODO: VT DEBUG
//        MatrixUtils.matrixMultiply(internalBofZ, preX,
//                globalThreadRowRange.getLength(), targetDimension,
//                ParallelOps.globalColCount, blockSize, outMM);
    }

    private void calculateBofZ(
            double[] preX, int targetDimension, double tCur, short[]
            distances, WeightsWrap1D weights,
            double[][] outBofZ) {

        int threadRowCount = globalThreadRowRange.getLength();

        double vBlockValue = -1;

        double diff = 0.0;
        if (tCur > 10E-10) {
            diff = Math.sqrt(2.0 * targetDimension) * tCur;
        }

        double[] outBofZLocalRow;
        double origD, weight, dist;

        final int globalColCount = config.numberDataPoints;
        final int globalRowOffset = globalThreadRowRange.getStartIndex();
        int globalRow;
        for (int threadLocalRow = 0; threadLocalRow < threadRowCount;
             ++threadLocalRow) {
            globalRow = threadLocalRow + globalRowOffset;
            outBofZLocalRow = outBofZ[threadLocalRow];
            outBofZLocalRow[globalRow] = 0;
            for (int globalCol = 0; globalCol < config.numberDataPoints;
                 globalCol++) {
                 /* B_ij = - w_ij * delta_ij / d_ij(Z), if (d_ij(Z) != 0) 0,
				 * otherwise v_ij = - w_ij.
				 *
				 * Therefore, B_ij = v_ij * delta_ij / d_ij(Z). 0 (if d_ij
				 * (Z) >=
				 * small threshold) --> the actual meaning is (if d_ij(Z) == 0)
				 * BofZ[i][j] = V[i][j] * deltaMat[i][j] / CalculateDistance
				 * (ref
				 * preX, i, j);*/

                // this is for the i!=j case. For i==j case will be calculated
                // separately (see above).
                if (globalRow == globalCol) continue;


                origD = distances[threadLocalRow * globalColCount +
                        globalCol] * INV_SHORT_MAX;
////                weight = weights.getWeight(threadLocalRow, globalCol);
                weight = 1.0;
//
                if (origD < 0 || weight == 0) {
                    continue;
                }

                dist = calculateEuclideanDist(preX, globalRow, globalCol,
                        targetDimension);
                if (dist >= 1.0E-10 && diff < origD) {
                    outBofZLocalRow[globalCol] = (weight * vBlockValue *
                            (origD - diff) / dist);
                } else {
                    outBofZLocalRow[globalCol] = 0;
                }
                outBofZLocalRow[globalRow] -= outBofZLocalRow[globalCol];
            }
        }
    }

    private static void extractPoints(
            Bytes bytes, int numPoints, int dimension, double[] to) {
        int pos = 0;
        int offset;
        for (int i = 0; i < numPoints; ++i) {
            offset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                bytes.position(pos);
                to[offset + j] = bytes.readDouble(pos);
                pos += Double.BYTES;
            }
        }
    }

    private static void mergePartials(double[][] partials, double[] result) {
        int offset = 0;
        for (double[] partial : partials) {
            System.arraycopy(partial, 0, result, offset, partial.length);
            offset += partial.length;
        }
    }

    private static void mergePartials(
            double[][] partials, Bytes result) {
        result.position(0);
        for (double[] partial : partials) {
            for (double aPartial : partial) {
                result.writeDouble(aPartial);
            }
        }
    }




    /*private static double calculateEuclideanDist(
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
    }*/

    public double calculateEuclideanDist(double[] v, int i, int j, int d) {
        double t = 0.0;
        double e;
        i = d * i;
        j = d * j;
        for (int k = 0; k < d; ++k) {
            e = v[i + k] - v[j + k];
            t += e * e;
        }
        return Math.sqrt(t);
    }





    private static TransformationFunction loadFunction(String classFile) {
        ClassLoader classLoader = ProgramLRT.class.getClassLoader();
        try {
            Class aClass = classLoader.loadClass(classFile);
            return (TransformationFunction) aClass.newInstance();
        } catch (ClassNotFoundException | InstantiationException |
                IllegalAccessException e) {
            throw new RuntimeException("Failed to load class: " + classFile, e);
        }
    }

    private void readDistancesAndWeights(boolean isSammon) {
        TransformationFunction function;
        if (!Strings.isNullOrEmpty(config.transformationFunction)) {
            function = loadFunction(config.transformationFunction);
        } else {
            function = (config.distanceTransform != 1.0
                    ? (d -> Math.pow(d, config.distanceTransform))
                    : null);
        }


        int elementCount = globalThreadRowRange.getLength() * config.numberDataPoints;
        distances = new short[elementCount];
        if (config.repetitions == 1) {
            BinaryReader1D.readRowRange(config.distanceMatrixFile,
                    globalThreadRowRange, config.numberDataPoints, byteOrder,
                    true, function, distances);
        } else {
            BinaryReader1D.readRowRange(config.distanceMatrixFile,
                    globalThreadRowRange, config.numberDataPoints, byteOrder,
                    true, function, config.repetitions, distances);
        }

        if (!Strings.isNullOrEmpty(config.weightMatrixFile)) {
            short[] w = null;
            function = !Strings.isNullOrEmpty(config
                    .weightTransformationFunction)
                    ? loadFunction(config.weightTransformationFunction)
                    : null;
            if (!config.isSimpleWeights) {
                w = new short[elementCount];
                if (config.repetitions == 1) {
                    BinaryReader1D.readRowRange(config.weightMatrixFile,
                            globalThreadRowRange, config.numberDataPoints,
                            byteOrder, true, function, w);
                } else {
                    BinaryReader1D.readRowRange(config.weightMatrixFile,
                            globalThreadRowRange, config.numberDataPoints,
                            byteOrder, true, function, config.repetitions, w);
                }
                weights = new WeightsWrap1D(
                        w, distances, isSammon, config.numberDataPoints);
            } else {
                double[] sw = null;
                sw = BinaryReader2D.readSimpleFile(config.weightMatrixFile,
                        config.numberDataPoints);
                weights = new WeightsWrap1D(sw, globalThreadRowRange,
                        distances, isSammon, config.numberDataPoints,
                        function);
            }
        } else {
            weights = new WeightsWrap1D(
                    null, distances, isSammon, config.numberDataPoints);
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
        } catch (ParseException e) {
            e.printStackTrace();
        }
        return Optional.fromNullable(null);
    }

}
