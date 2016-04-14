package edu.indiana.soic.spidal.damds;

import com.google.common.base.Optional;
import com.google.common.base.Strings;
import edu.indiana.soic.spidal.common.*;
import edu.indiana.soic.spidal.configuration.section.DAMDSSection;
import edu.indiana.soic.spidal.damds.threads.ThreadCommunicator;
import edu.indiana.soic.spidal.damds.timing.*;
import mpi.MPI;
import mpi.MPIException;
import net.openhft.lang.io.Bytes;
import org.apache.commons.cli.*;

import java.io.*;
import java.nio.ByteOrder;
import java.nio.LongBuffer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.BrokenBarrierException;
import java.util.regex.Pattern;
import java.util.stream.IntStream;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;

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

    private double[] partialSigma;
    private double[] vArray;

    //Config Settings
    private DAMDSSection config;
    private ByteOrder byteOrder;
    private short[] distances;
    private WeightsWrap1D weights;

    private int BlockSize;

    private int threadId;
    private Range threadRowRange;

    final private RefObj<Integer> refInt = new RefObj<>();

    private ThreadCommunicator comm;

    public ProgramWorker(int threadId, ThreadCommunicator comm, DAMDSSection config, ByteOrder byteOrder, int blockSize){
        System.out.println("In ProgramWorker Constructor " + threadId);
        this.threadId = threadId;
        this.comm = comm;
        this.config = config;
        this.byteOrder = byteOrder;
        this.BlockSize = blockSize;
    }

    public void setup() {
        threadRowRange = ParallelOps.threadRowRanges[threadId];
    }
    public void  run() {
        try {
            setup();

            readDistancesAndWeights(config.isSammon);

            RefObj<Integer> missingDistCount = new RefObj<>();
            DoubleStatistics distanceSummary = calculateStatistics(
                distances, weights, missingDistCount);
            System.out.println("*****CAME HERE****");
            double missingDistPercent = missingDistCount.getValue() /
                                        (Math.pow(config.numberDataPoints, 2));
            INV_SUM_OF_SQUARE = 1.0/distanceSummary.getSumOfSquare();
            Utils.printMessage(
                "\nDistance summary... \n" + distanceSummary.toString() +
                "\n  MissingDistPercentage=" +
                missingDistPercent);

            weights.setAvgDistForSammon(distanceSummary.getAverage());
            changeZeroDistancesToPostiveMin(distances, distanceSummary.getPositiveMin());

            /*
            // Allocating point arrays once for all
            allocateArrays();

            if (Strings.isNullOrEmpty(config.initialPointsFile)){
                generateInitMapping(
                    config.numberDataPoints, config.targetDimension, preX);
            } else {
                readInitMapping(config.initialPointsFile, preX, config.targetDimension);
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

                *//* Note - quick way to test programs without running full
                * number of temperature loops *//*
                if (config.maxtemploops > 0 && loopNum == config.maxtemploops){
                    break;
                }
            }
            loopTimer.stop();*/
        }
        catch (MPIException e) {
            Utils.printAndThrowRuntimeException(new RuntimeException(e));
        }
        catch (InterruptedException e) {
            e.printStackTrace();
        }
        catch (BrokenBarrierException e) {
            e.printStackTrace();
        }
    }

    /*
    private static void allocateArrays() {
        // Allocating point arrays once for all
        final int numberDataPoints = config.numberDataPoints;
        final int targetDimension = config.targetDimension;
        // Note - prex[][] to prex[]
        //        preX = new double[numberDataPoints][targetDimension];
        preX = new double[numberDataPoints*targetDimension];

        // Note - BC[][] to BC[]
        //        BC = new double[numberDataPoints][targetDimension];
        BC = new double[numberDataPoints*targetDimension];

        // Note - MMr[][] to MMr[]
        MMr = new double[numberDataPoints * targetDimension];
        // Note - MMAp[][] to MMAp[]
        //        MMAp = new double[numberDataPoints][targetDimension];
        MMAp = new double[numberDataPoints*targetDimension];
        final int threadCount = ParallelOps.threadCount;
        threadPartialBofZ = new double[threadCount][][];
        // Note - threadPartialMM[][][] to threadPartialMM[][]
        //        threadPartialMM = new double[threadCount][][];
        threadPartialMM = new double[threadCount][];
        vArray = new double[threadCount][];
        int threadRowCount;
        for (int i = 0; i < threadCount; ++i){
            threadRowCount = ParallelOps.threadRowCounts[i];
            threadPartialBofZ[i] = new double[threadRowCount][ParallelOps.globalColCount];
            threadPartialMM[i] = new double[threadRowCount * config.targetDimension];
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

    */

    private void changeZeroDistancesToPostiveMin(
        short[] distances, double positiveMin) {
        double tmpD;
        for (int i = 0; i < distances.length; ++i){
            tmpD = distances[i] * INV_SHORT_MAX;
            if (tmpD < positiveMin && tmpD >= 0.0){
                distances[i] = (short)(positiveMin * SHORT_MAX);
            }
        }
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
        String initialPointsFile, double[] preX, int dimension) {
        try (BufferedReader br = Files.newBufferedReader(Paths.get(initialPointsFile),
            Charset
                .defaultCharset())){
            String line;
            Pattern pattern = Pattern.compile("[\t]");
            int row = 0;
            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line))
                    continue; // continue on empty lines - "while" will break on null anyway;

                String[] splits = pattern.split(line.trim());

                for (int i = 0; i < dimension; ++i){
                    preX[row+i] = Double.parseDouble(splits[i].trim());
                }
                row+=dimension;
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
            /*writer.println(labels.get(String.valueOf(ids[i]))); // print label*/
            // value, which
            // is
            // ONE for all data.
        }
        writer.flush();
        writer.close();
    }

    /*private void generateVArray(
        short[] distances, WeightsWrap1D weights, double[][] vArray) {

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
    }*/

    /*private static void generateVArrayInternal(
        Integer threadIdx, short[] distances, WeightsWrap1D weights, double[] v) {
        int threadRowCount = ParallelOps.threadRowCounts[threadIdx];

        int rowOffset = ParallelOps.threadRowStartOffsets[threadIdx] +
                        ParallelOps.procRowStartOffset;
        for (int i = 0; i < threadRowCount; ++i) {
            int globalRow = i + rowOffset;
            int procLocalRow = globalRow - ParallelOps.procRowStartOffset;
            for (int globalCol = 0; globalCol < ParallelOps.globalColCount; ++globalCol) {
                if (globalRow == globalCol) continue;

                double origD = distances[procLocalRow*ParallelOps.globalColCount+globalCol] * INV_SHORT_MAX;
                double weight = weights.getWeight(procLocalRow, globalCol);

                if (origD < 0 || weight == 0) {
                    continue;
                }

                v[i] += weight;
            }
            v[i] += 1;
        }
    }*/

    /*private static void calculateConjugateGradient(
        double[] preX, int targetDimension, int numPoints, double[] BC, int cgIter, double cgThreshold,
        RefObj<Integer> outCgCount, RefObj<Integer> outRealCGIterations,
        WeightsWrap1D weights, int blockSize, double[][] vArray, double[] MMr, double[] MMAp, double[][] threadPartialMM)

        throws MPIException {


        zeroOutArray(threadPartialMM);
        CGTimings.startTiming(CGTimings.TimingTask.MM);
        calculateMM(preX, targetDimension, numPoints, weights, blockSize,
            vArray, MMr, threadPartialMM);

        CGTimings.endTiming(CGTimings.TimingTask.MM);
        // This barrier was necessary for correctness when using
        // a single mmap file
        ParallelOps.worldProcsComm.barrier();

        int iOffset;
        double[] tmpRHSRow;

        for(int i = 0; i < numPoints; ++i) {
            iOffset = i*targetDimension;
            for (int j = 0; j < targetDimension; ++j) {
                BC[iOffset+j] -= MMr[iOffset+j];
                MMr[iOffset+j] = BC[iOffset+j];
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
                iOffset = i*targetDimension;
                for (int j = 0; j < targetDimension; ++j) {
                    preX[iOffset+j] += alpha * BC[iOffset+j];
                }
            }

            if (rTr < testEnd) {
                break;
            }

            //update ri to ri+1
            for(int i = 0; i < numPoints; ++i) {
                iOffset = i*targetDimension;
                for (int j = 0; j < targetDimension; ++j) {
                    MMr[iOffset+j] -= alpha * MMAp[iOffset+j];
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
                iOffset = i*targetDimension;
                for (int j = 0; j < targetDimension; ++j) {
                    BC[iOffset+j] = MMr[iOffset+j] + beta * BC[iOffset+j];
                }
            }

        }
        CGTimings.endTiming(CGTimings.TimingTask.CG_LOOP);
        outCgCount.setValue(outCgCount.getValue() + cgCount);
    }*/

    private static void calculateMM(
        double[] x, int targetDimension, int numPoints, WeightsWrap1D weights,
        int blockSize, double[][] vArray, double[] outMM,
        double[][] internalPartialMM) throws MPIException {

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
            mergePartials(internalPartialMM, ParallelOps.mmapXWriteBytes);
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
            mergePartials(internalPartialMM, outMM);
        }
    }

    private static void calculateMMInternal(
        Integer threadIdx, double[] x, int targetDimension, int numPoints,
        WeightsWrap1D weights, int blockSize, double[][] vArray, double[] outMM) {

        MatrixUtils
            .matrixMultiplyWithThreadOffset(weights, vArray[threadIdx], x,
                ParallelOps.threadRowCounts[threadIdx], targetDimension,
                numPoints, blockSize,
                ParallelOps.threadRowStartOffsets[threadIdx],
                ParallelOps.threadRowStartOffsets[threadIdx]
                + ParallelOps.procRowStartOffset, outMM);
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

    /*private static void calculateBC(
        double[] preX, int targetDimension, double tCur, short[] distances,
        WeightsWrap1D weights, int blockSize, double[] BC,
        double[][][] threadPartialBCInternalBofZ,
        double[][] threadPartialBCInternalMM)
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
                0, preX, targetDimension, tCur, distances, weights, blockSize,
                threadPartialBCInternalBofZ[0], threadPartialBCInternalMM[0]);
            BCTimings.endTiming(
                BCTimings.TimingTask.BC_INTERNAL, 0);
        }

        if (ParallelOps.worldProcsCount > 1) {
            BCTimings.startTiming(BCTimings.TimingTask.BC_MERGE, 0);
            mergePartials(threadPartialBCInternalMM, ParallelOps.mmapXWriteBytes);
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
            mergePartials(threadPartialBCInternalMM, BC);
        }
    }*/

    /*private static void calculateBCInternal(
        Integer threadIdx, double[] preX, int targetDimension, double tCur,
        short[] distances, WeightsWrap1D weights, int blockSize,
        double[][] internalBofZ, double[] outMM) {

        BCInternalTimings.startTiming(BCInternalTimings.TimingTask.BOFZ, threadIdx);
        calculateBofZ(threadIdx, preX, targetDimension, tCur,
            distances, weights, internalBofZ);
        BCInternalTimings.endTiming(BCInternalTimings.TimingTask.BOFZ, threadIdx);

        // Next we can calculate the BofZ * preX.
        BCInternalTimings.startTiming(BCInternalTimings.TimingTask.MM, threadIdx);
        MatrixUtils.matrixMultiply(internalBofZ, preX,
            ParallelOps.threadRowCounts[threadIdx], targetDimension,
            ParallelOps.globalColCount, blockSize, outMM);
        BCInternalTimings.endTiming(BCInternalTimings.TimingTask.MM, threadIdx);
    }*/

    /*private static void calculateBofZ(
        int threadIdx, double[] preX, int targetDimension, double tCur, short[] distances, WeightsWrap1D weights,
        double[][] outBofZ) {

        int threadRowCount = ParallelOps.threadRowCounts[threadIdx];

        double vBlockValue = -1;

        double diff = 0.0;
        if (tCur > 10E-10) {
            diff = Math.sqrt(2.0 * targetDimension)  * tCur;
        }

        short[] distancesProcLocalRow;
        double[] outBofZLocalRow;
        double[] preXGlobalRow;
        double origD, weight, dist;

        final int globalColCount = ParallelOps.globalColCount;
        final int globalRowOffset = ParallelOps.threadRowStartOffsets[threadIdx]
                                    + ParallelOps.procRowStartOffset;
        int globalRow, procLocalRow;
        for (int localRow = 0; localRow < threadRowCount; ++localRow) {
            globalRow = localRow + globalRowOffset;
            procLocalRow = globalRow - ParallelOps.procRowStartOffset;
            outBofZLocalRow = outBofZ[localRow];
            outBofZLocalRow[globalRow] = 0;
            for (int globalCol = 0; globalCol < ParallelOps.globalColCount; globalCol++) {
				*//*
				 * B_ij = - w_ij * delta_ij / d_ij(Z), if (d_ij(Z) != 0) 0,
				 * otherwise v_ij = - w_ij.
				 *
				 * Therefore, B_ij = v_ij * delta_ij / d_ij(Z). 0 (if d_ij(Z) >=
				 * small threshold) --> the actual meaning is (if d_ij(Z) == 0)
				 * BofZ[i][j] = V[i][j] * deltaMat[i][j] / CalculateDistance(ref
				 * preX, i, j);
				 *//*
                // this is for the i!=j case. For i==j case will be calculated
                // separately (see above).
                if (globalRow == globalCol) continue;


                origD = distances[procLocalRow*globalColCount+globalCol] * INV_SHORT_MAX;
                weight = weights.getWeight(procLocalRow,globalCol);

                if (origD < 0 || weight == 0) {
                    continue;
                }

                dist = calculateEuclideanDist(preX, globalRow, globalCol,targetDimension);
                if (dist >= 1.0E-10 && diff < origD) {
                    outBofZLocalRow[globalCol] = (weight * vBlockValue * (origD - diff) / dist);
                } else {
                    outBofZLocalRow[globalCol] = 0;
                }

                outBofZLocalRow[globalRow] -= outBofZLocalRow[globalCol];
            }
        }
    }*/

    private static void extractPoints(
        Bytes bytes, int numPoints, int dimension, double[] to) {
        int pos = 0;
        int offset;
        for (int i = 0; i < numPoints; ++i){
            offset = i*dimension;
            for (int j = 0; j < dimension; ++j) {
                bytes.position(pos);
                to[offset+j] = bytes.readDouble(pos);
                pos += Double.BYTES;
            }
        }
    }

    private static void mergePartials(double[][] partials, double[] result){
        int offset = 0;
        for (double [] partial : partials){
            System.arraycopy(partial, 0, result, offset, partial.length);
            offset+=partial.length;
        }
    }

    private static void mergePartials(
        double[][] partials, Bytes result){
        result.position(0);
        for (double [] partial : partials){
            for (double aPartial : partial) {
                result.writeDouble(aPartial);
            }
        }
    }

    private double calculateStress(
        double[] preX, int targetDimension, double tCur, short[] distances,
        WeightsWrap1D weights, double invSumOfSquareDist, double[] internalPartialSigma)
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
            StressTimings.startTiming(StressTimings.TimingTask.STRESS_INTERNAL, 0);
            stress = calculateStressInternal(0, preX, targetDimension, tCur,
                distances, weights);
            StressTimings.endTiming(StressTimings.TimingTask.STRESS_INTERNAL, 0);
        }

        if (ParallelOps.worldProcsCount > 1) {
            // Write thread local reduction to shared memory map
            ParallelOps.mmapSWriteBytes.position(0);
            ParallelOps.mmapSWriteBytes.writeDouble(stress);

            // Important barrier here - as we need to make sure writes are done
            // to the mmap file.
            // It's sufficient to wait on ParallelOps.mmapProcComm,
            // but it's cleaner for timings if we wait on the whole world
            ParallelOps.worldProcsComm.barrier();
            if (ParallelOps.isMmapLead) {
                // Node local reduction using shared memory maps
                ParallelOps.mmapSReadBytes.position(0);
                stress = 0.0;
                for (int i = 0; i < ParallelOps.mmapProcsCount; ++i){
                    stress += ParallelOps.mmapSReadBytes.readDouble();
                }
                ParallelOps.mmapSWriteBytes.position(0);
                ParallelOps.mmapSWriteBytes.writeDouble(stress);

                // Leaders participate in MPI AllReduce
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
        }
        return stress * invSumOfSquareDist;
    }

    private double calculateStressInternal(
        int threadIdx, double[] preX, int targetDim, double tCur, short[] distances, WeightsWrap1D weights) {

        StressInternalTimings.startTiming(StressInternalTimings.TimingTask.COMP, threadIdx);
        double sigma = 0.0;
        double diff = 0.0;
        if (tCur > 10E-10) {
            diff = Math.sqrt(2.0 * targetDim) * tCur;
        }

        int threadRowCount = ParallelOps.threadRowCounts[threadIdx];
        final int globalRowOffset = ParallelOps.threadRowStartOffsets[threadIdx]
                                    + ParallelOps.procRowStartOffset;

        int globalColCount = ParallelOps.globalColCount;
        int globalRow, procLocalRow;
        double origD, weight, euclideanD;
        double heatD, tmpD;
        for (int localRow = 0; localRow < threadRowCount; ++localRow){
            globalRow = localRow + globalRowOffset;
            procLocalRow = globalRow - ParallelOps.procRowStartOffset;
            for (int globalCol = 0; globalCol < globalColCount; globalCol++) {
                origD = distances[procLocalRow * globalColCount + globalCol]
                        * INV_SHORT_MAX;
                weight = weights.getWeight(procLocalRow,globalCol);

                if (origD < 0 || weight == 0) {
                    continue;
                }

                euclideanD = globalRow != globalCol ? calculateEuclideanDist(
                    preX, globalRow , globalCol, targetDim) : 0.0;

                heatD = origD - diff;
                tmpD = origD >= diff ? heatD - euclideanD : -euclideanD;
                sigma += weight * tmpD * tmpD;
            }
        }
        StressInternalTimings.endTiming(StressInternalTimings.TimingTask.COMP, threadIdx);
        return sigma;
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

    public static double calculateEuclideanDist(double[] v, int i, int j, int d){
        double t = 0.0; double e;
        i = d*i; j=d*j;
        for (int k = 0; k < d; ++k){
            e = v[i+k] - v[j+k];
            t += e*e;
        }
        return Math.sqrt(t);
    }

    private static void generateInitMapping(
        int numPoints, int targetDim, double[] preX) throws MPIException {

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

    private DoubleStatistics calculateStatistics(
        short[] distances, WeightsWrap1D weights, RefObj<Integer> missingDistCount)
        throws MPIException, BrokenBarrierException, InterruptedException {

        DoubleStatistics distanceSummary =
            calculateStatisticsInternal(distances, weights, refInt);
        /*comm.sumOverThreads(threadId, distanceSummary);
        comm.sumOverThreads(threadId, refInt);*/

        if (ParallelOps.worldProcsCount > 1 && threadId == 0) {
            distanceSummary = ParallelOps.allReduce(distanceSummary);
            refInt.setValue(ParallelOps.allReduce(refInt.getValue()));
        }
        /*comm.threadBarrier();
        comm.bcastOverThreads(threadId, distanceSummary, 0);
        comm.bcastOverThreads(threadId, refInt, 0);*/
        missingDistCount.setValue(refInt.getValue());
        return distanceSummary;
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

    private void readDistancesAndWeights(boolean isSammon) {
        TransformationFunction function;
        if (!Strings.isNullOrEmpty(config.transformationFunction)) {
            function = loadFunction(config.transformationFunction);
        } else {
            function = (config.distanceTransform != 1.0
                            ? (d -> Math.pow(d, config.distanceTransform))
                            : null);
        }


        int elementCount = threadRowRange.getLength() * ParallelOps.globalColCount;
        distances = new short[elementCount];
        if (config.repetitions == 1){
            BinaryReader1D.readRowRange(config.distanceMatrixFile,
                threadRowRange, ParallelOps.globalColCount, byteOrder,
                true, function, distances);
        }else{
            BinaryReader1D.readRowRange(config.distanceMatrixFile,
                threadRowRange, ParallelOps.globalColCount, byteOrder,
                true, function, config.repetitions, distances);
        }

        if (!Strings.isNullOrEmpty(config.weightMatrixFile)){
            short[] w = null;
            function = !Strings.isNullOrEmpty(config.weightTransformationFunction)
                ? loadFunction(config.weightTransformationFunction)
                : null;
            if (!config.isSimpleWeights) {
                w = new short[elementCount];
                if (config.repetitions == 1) {
                    BinaryReader1D.readRowRange(config.weightMatrixFile,
                        threadRowRange, ParallelOps.globalColCount,
                        byteOrder, true, function, w);
                }
                else {
                    BinaryReader1D.readRowRange(config.weightMatrixFile,
                        threadRowRange, ParallelOps.globalColCount,
                        byteOrder, true, function, config.repetitions, w);
                }
                weights = new WeightsWrap1D(
                    w, distances, isSammon, ParallelOps.globalColCount);
            } else {
                double[] sw = null;
                sw = BinaryReader2D.readSimpleFile(config.weightMatrixFile, config.numberDataPoints);
                weights = new WeightsWrap1D(sw, threadRowRange, distances, isSammon, ParallelOps.globalColCount, function);
            }
        } else {
            weights = new WeightsWrap1D(
                null, distances, isSammon, ParallelOps.globalColCount);
        }

    }

    private DoubleStatistics calculateStatisticsInternal(
        short[] distances, WeightsWrap1D weights, RefObj<Integer> refMissingDistCount) {

        int missingDistCount = 0;
        DoubleStatistics stat = new DoubleStatistics();

        int threadRowCount = ParallelOps.threadRowCounts[threadId];

        double origD, weight;
        System.out.println("Thread Row Start Offset " + ParallelOps.threadRowStartOffsets[threadId]);
        for (int localRow = 0; localRow < threadRowCount; ++localRow){
            for (int globalCol = 0; globalCol < ParallelOps.globalColCount; globalCol++) {
                origD = distances[localRow*ParallelOps.globalColCount + globalCol] * INV_SHORT_MAX;
                weight = weights.getWeight(localRow,globalCol);
                if (origD < 0) {
                    // Missing distance
                    ++missingDistCount;
                    continue;
                }
                if (weight == 0) continue; // Ignore zero weights

                stat.accept(origD);
            }
        }
        refMissingDistCount.setValue(missingDistCount);
        System.out.println("In calcstatinternal " + threadId + " threadRowCount=" + threadRowCount);
        return stat;
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
