package edu.indiana.soic.spidal.damds;

import com.google.common.base.Optional;
import com.google.common.base.Stopwatch;
import com.google.common.base.Strings;
import edu.indiana.soic.spidal.common.BinaryReader;
import edu.indiana.soic.spidal.common.DoubleStatistics;
import edu.indiana.soic.spidal.common.MatrixUtils;
import edu.indiana.soic.spidal.common.RefObj;
import edu.indiana.soic.spidal.configuration.ConfigurationMgr;
import edu.indiana.soic.spidal.configuration.section.DAMDSSection;
import mpi.MPIException;
import org.apache.commons.cli.*;

import java.io.*;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
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
    }

    //Config Settings
    public static DAMDSSection config;
    public static ByteOrder byteOrder;
    public static BinaryReader distances;
    public static BinaryReader weights;

    public static int BlockSize = 64;
    public static double [][] distanceArray;
    public static double [][] weightsArray;

    /* TODO fix these variables to probably local variables and use ref passing*/
    public static int CG_REAL_ITER;
    public static int SMACOF_REAL_ITER;
    /* This is the maximum determinant value for using Inverse or pseudo inverse*/
    public static double MAX_DET = 10000000;
    public static int MAX_ITER = 10000;

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
    public static void main(String[] args) {
        Stopwatch mainTimer = Stopwatch.createStarted();
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

            Utils.printMessage("\n== DAMDS run started on " + new Date() + " ==\n");
            Utils.printMessage(config.toString(false));

            readDistancesAndWeights();
            RefObj<Integer> missingDistCount = new RefObj<>();
            DoubleStatistics distanceSummary = calculateStatistics(
                missingDistCount);

            double missingDistPercent = missingDistCount.getValue() /
                                        (Math.pow(config.numberDataPoints, 2));
            Utils.printMessage(
                "\nDistance summary... \n" + distanceSummary.toString() +
                "\n  MissingDistPercentage=" +
                missingDistPercent);

            double [][] vArray = generateVArray();

            double[][] preX = generateInitMapping(
                config.numberDataPoints, config.targetDimension);
            double tCur = 0.0;
            double preStress = calculateStress(
                preX, tCur, config.targetDimension, config.isSammon,
                distanceSummary.getAverage(), distanceSummary.getSumOfSquare());
            Utils.printMessage("\nInitial stress=" + preStress);

            double X[][] = null;
            double BC[][] = null;

            double tMax =
                calculateMaxT(distanceSummary.getMax(), config.targetDimension);
            double tMin = (0.01 * tMax < 0.01) ? 0.01 * tMax : 0.01;
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
            while (true) {
                preStress = calculateStress(
                    preX, tCur, config.targetDimension, config.isSammon,
                    distanceSummary.getAverage(),
                    distanceSummary.getSumOfSquare());
                diffStress = config.threshold + 1.0;

                Utils.printMessage(
                    String.format(
                        "\nStart of loop %d Temperature (T_Cur) %.5g",
                        loopNum, tCur));

                int itrNum = 0;
                RefObj<Integer> cgCount = new RefObj<>(0);
                while (diffStress >= config.threshold) {
                    BC = calculateBC(
                        preX, config.targetDimension,
                        tCur, config.isSammon, distanceSummary.getAverage(),
                        BlockSize);

                    X = calculateConjugateGradient(
                        preX, config.targetDimension, config.numberDataPoints,
                        BC, config.cgIter, config.cgErrorThreshold, cgCount,
                        config.isSammon, distanceSummary.getAverage(),
                        BlockSize, vArray);

                    stress = calculateStress(
                        X, tCur, config.targetDimension, config.isSammon,
                        distanceSummary.getAverage(),
                        distanceSummary.getSumOfSquare());

                    diffStress = preStress - stress;
                    preStress = stress;
                    preX = MatrixUtils.copy(X);

                    if ((itrNum % 10 == 0) || (itrNum >= MAX_ITER)) {
                        Utils.printMessage(
                            String.format(
                                "  Loop %d Iteration %d Avg CG count %.5g " +
                                "Stress " +
                                "%.5g", loopNum, itrNum,
                                (cgCount.getValue() * 1.0 / (itrNum + 1)),
                                stress));
                    }
                    ++itrNum;
                    ++SMACOF_REAL_ITER;
                }
                --itrNum;
                if (itrNum >=0 && !(itrNum % 10 == 0) && !(itrNum >=
                                                           MAX_ITER)) {
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


            /* TODO Fix error handling here */
            if (Strings.isNullOrEmpty(config.labelFile) || config.labelFile.toUpperCase().endsWith(
                "NOLABEL")) {
                try {
                    writeOuput(X, config.pointsFile);
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
            } else {
                try {
                    writeOuput(X, config.labelFile, config.pointsFile);
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
            }

            Double finalStress = calculateStress(X,tCur, config.targetDimension, config.isSammon,
                                                 distanceSummary.getAverage(),
                                                 distanceSummary.getSumOfSquare());

            mainTimer.stop();


            Utils.printMessage("Finishing DAMDS run ...");
            long totalTime = mainTimer.elapsed(TimeUnit.MILLISECONDS);
            long loopTime = loopTimer.elapsed(TimeUnit.MILLISECONDS);
            Utils.printMessage(
                String.format(
                    "  Total Time: %s (%d ms) Loop Time: %s (%d ms)",
                    formatElapsedMillis(totalTime), totalTime,
                    formatElapsedMillis(loopTime), loopTime));
            Utils.printMessage("  Total Loops: " + loopNum);
            Utils.printMessage("  Total Iterations: " + SMACOF_REAL_ITER);
            Utils.printMessage(
                String.format(
                    "  Total CG Iterations: %d Avg. CG Iterations: %.5g",
                    CG_REAL_ITER, (CG_REAL_ITER * 1.0) / SMACOF_REAL_ITER));
            Utils.printMessage("  Final Stress:\t" + finalStress);
            Utils.printMessage("== DAMDS run completed on " + new Date() + " ==");

            ParallelOps.tearDownParallelism();
        }
        catch (MPIException e) {
            Utils.printAndThrowRuntimeException(new RuntimeException(e));
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
            writer.print(String.valueOf(i) + "\t"); // print ID.
            for (int j = 0; j < vecLen; j++) {
                writer.print(format.format(x[i][j]) + "\t"); // print
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
        String line = null;
        String parts[] = null;
        Map<String, Integer> labels = new HashMap<String, Integer>();
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
            writer.print(String.valueOf(i) + "\t"); // print ID.
            for (int j = 0; j < vecLen; j++) {
                writer.print(format.format(X[i][j]) + "\t"); // print
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

    private static double[][] generateVArray() {
        double [][] vArray = new double[ParallelOps.threadCount][];

        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> vArray[threadIdx] =
                        generateVArrayInternal(threadIdx)));
        }
        else {
            vArray[0] = generateVArrayInternal(0);
        }
        return vArray;
    }

    private static double[] generateVArrayInternal(Integer threadIdx) {
        int threadRowCount = ParallelOps.threadRowCounts[threadIdx];
        double [] v = new double[threadRowCount];

        int rowOffset = ParallelOps.threadRowStartOffsets[threadIdx] +
                     ParallelOps.procRowStartOffset;
        int threadPointStartOffset =
            ParallelOps.threadPointStartOffsets[threadIdx];
        for (int i = 0; i < threadRowCount; ++i) {
            int globalRow = i + rowOffset;
            for (int j = 0; j < ParallelOps.globalColCount; ++j) {
                if (globalRow == j) continue;

                int procLocalPnum = ((i * ParallelOps.globalColCount) + j) + threadPointStartOffset;
                double origD = distances.getValue(procLocalPnum);
                double weight = weights.getValue(procLocalPnum);

                if (origD < 0 || weight == 0) {
                    continue;
                }

                v[i] += weight;
            }
            v[i] += 1;
        }
        return v;
    }

    private static double[][] calculateConjugateGradient(
        double[][] preX, int targetDimension, int numPoints, double[][] BC, int cgIter, double cgThreshold,
        RefObj<Integer> outCgCount, boolean isSammon, double avgDist, int blockSize, double [][] vArray)

        throws MPIException {

        double[][] X = null;
        double[][] r = new double[numPoints][targetDimension];
        double[][] p = new double[numPoints][targetDimension];

        X = preX;
        r = calculateMM(X, targetDimension, numPoints, isSammon, avgDist, blockSize,
                        vArray);

        for(int i = 0; i < numPoints; ++i)
            for(int j = 0; j < targetDimension; ++j){
                p[i][j] = BC[i][j] - r[i][j];
                r[i][j] = p[i][j];
            }

        int cgCount = 0;
        double rTr = innerProductCalculation(r);
        // Adding relative value test for termination as suggested by Dr. Fox.
        double testEnd = rTr * cgThreshold;

        //System.out.println("1");
        while(cgCount < cgIter){
            cgCount++;
            ++CG_REAL_ITER;
            //System.out.println("2");
            //calculate alpha
            double[][] Ap = calculateMM(p, targetDimension, numPoints, isSammon, avgDist, blockSize, vArray);

            double alpha = rTr
                           /innerProductCalculation(p, Ap);

            //update Xi to Xi+1
            for(int i = 0; i < numPoints; ++i)
                for(int j = 0; j < targetDimension; ++j)
                    X[i][j] += alpha * p[i][j];

            if (rTr < testEnd) {
                break;
            }

            //update ri to ri+1
            for(int i = 0; i < numPoints; ++i)
                for(int j = 0; j < targetDimension; ++j)
                    r[i][j] = r[i][j] - alpha * Ap[i][j];

            //calculate beta
            double rTr1 = innerProductCalculation(r);
            double beta = rTr1/rTr;
            rTr = rTr1;

            //update pi to pi+1
            for(int i = 0; i < numPoints; ++i)
                for(int j = 0; j < targetDimension; ++j)
                    p[i][j] = r[i][j] + beta * p[i][j];


        }
        outCgCount.setValue(outCgCount.getValue() + cgCount);
        //		System.out.println("CGCount: " + cgCount + " TestEnd: " + testEnd + " rTr: " + rTr);
        return X;

    }

    private static double[][] calculateMM(
        double[][] x, int targetDimension, int numPoints, boolean isSammon, double avgDist,
        int blockSize, double[][] vArray) throws MPIException {

        double [][][] partialMMs = new double[ParallelOps.threadCount][][];

        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> partialMMs[threadIdx] =
                        calculateMMInternal(
                            threadIdx, x, targetDimension, numPoints, isSammon, avgDist, blockSize, vArray)));
        }
        else {
            partialMMs[0] = calculateMMInternal(
                0, x, targetDimension, numPoints, isSammon, avgDist, blockSize,
                vArray);
        }

        if (ParallelOps.procCount > 1) {
            mergePartials(partialMMs, targetDimension, ParallelOps.partialPointBuffer);
            return extractPoints(
                ParallelOps.allGather(
                    ParallelOps.partialPointBuffer, targetDimension),
                ParallelOps.globalColCount, targetDimension);
        } else {
            double [][] mm = new double[ParallelOps.globalColCount][targetDimension];
            mergePartials(partialMMs, targetDimension, mm);
            return mm;
        }
    }

    private static double[][] calculateMMInternal(
        Integer threadIdx, double[][] x, int targetDimension, int numPoints,
        boolean isSammon, double avgDist, int blockSize, double[][] vArray) {

        return MatrixUtils.matrixMultiply(
            (threadLocalRow, globalCol) -> {
                int procLocalPnum =
                    (threadLocalRow * ParallelOps.globalColCount) + globalCol +
                    ParallelOps.threadPointStartOffsets[threadIdx];
                double distance = distances.getValue(procLocalPnum);
                double weight = weights.getValue(procLocalPnum);
                return (distance < 0 || weight == 0) ? 0.0 :
                       (isSammon ? (1.0 /Math.max(distance,0.001 *avgDist)) : weight);
            }, vArray[threadIdx], x, ParallelOps.threadRowCounts[threadIdx],
            targetDimension, numPoints, blockSize,
            ParallelOps.threadRowStartOffsets[threadIdx] +
            ParallelOps.procRowStartOffset);
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

    private static double[][] calculateBC(
        double[][] preX, int targetDimension, double tCur, boolean isSammon,
        double avgDist, int blockSize) throws MPIException{

        double [][][] partialBCs = new double[ParallelOps.threadCount][][];

        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> partialBCs[threadIdx] =
                        calculateBCInternal(
                            threadIdx, preX, targetDimension, tCur, isSammon, avgDist, blockSize)));
        }
        else {
            partialBCs[0] = calculateBCInternal(
                0, preX, targetDimension, tCur, isSammon, avgDist, blockSize);
        }

        if (ParallelOps.procCount > 1) {
            mergePartials(partialBCs, targetDimension, ParallelOps.partialPointBuffer);
            return extractPoints(
                ParallelOps.allGather(
                    ParallelOps.partialPointBuffer, targetDimension),
                ParallelOps.globalColCount, targetDimension);
        } else {
            double [][] bc = new double[ParallelOps.globalColCount][targetDimension];
            mergePartials(partialBCs, targetDimension, bc);
            return bc;
        }
    }

    private static double[][] calculateBCInternal(
        Integer threadIdx, double[][] preX, int targetDimension, double tCur,
        boolean isSammon, double avgDist, int blockSize) {

        float [][] BofZ = calculateBofZ(threadIdx, preX, targetDimension, tCur, isSammon, avgDist);

        // Next we can calculate the BofZ * preX.
        return MatrixUtils.matrixMultiply(BofZ, preX, ParallelOps.threadRowCounts[threadIdx],
                                                  targetDimension, ParallelOps.globalColCount, blockSize);
    }

    private static float[][] calculateBofZ(int threadIdx, double[][] preX, int targetDimension, double tCur, boolean isSammon, double avgDist) {

        int threadRowCount = ParallelOps.threadRowCounts[threadIdx];
        float [][] BofZ = new float[threadRowCount][ParallelOps.globalColCount];

        double vBlockValue = (double) -1;

        double diff = 0.0;
        if (tCur > 10E-10) {
            diff = Math.sqrt(2.0 * targetDimension)  * tCur;
        }

        for (int localRow = 0; localRow < threadRowCount; ++localRow) {
            int globalRow = localRow + ParallelOps.threadRowStartOffsets[threadIdx] +
                     ParallelOps.procRowStartOffset;
            BofZ[localRow][globalRow] = 0;
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

                int procLocalPnum = (localRow + ParallelOps.threadRowStartOffsets[threadIdx]) * ParallelOps.globalColCount + globalCol;
                double origD = distances.getValue(procLocalPnum);
                double weight = weights.getValue(procLocalPnum);

                if (origD < 0 || weight == 0) {
                    continue;
                }

                weight = isSammon ? 1.0 / Math.max(origD, 0.001 * avgDist) : weight;

                double dist = calculateEuclideanDist(
                    preX, targetDimension, globalRow, globalCol);
                if (dist >= 1.0E-10 && diff < origD) {
                    BofZ[localRow][globalCol] = (float) (weight * vBlockValue * (origD - diff) / dist);
                } else {
                    BofZ[localRow][globalCol] = 0;
                }

                BofZ[localRow][globalRow] += -BofZ[localRow][globalCol];
            }
        }
        return BofZ;
    }

    private static double[][] extractPoints(
        DoubleBuffer buffer, int numPoints, int dimension) {
        double [][] points = new double[numPoints][dimension];
        int pos = 0;
        for (int i = 0; i < numPoints; ++i){
            buffer.position(pos);
            buffer.get(points[i]);
            pos += dimension;
        }
        return  points;
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

    private static void mergePartials(double [][][] partials, int dimension, DoubleBuffer result){
        int pos = 0;
        for (double [][] partial : partials){
            for (double [] point : partial){
                result.position(pos);
                result.put(point);
                pos += dimension;
            }
        }
    }

    private static double calculateMaxT(double maxOrigDistance, int targetDim) {
        double divider = Math.sqrt(2.0 * targetDim);
        return maxOrigDistance / divider;
    }

    private static double calculateStress(double[][] preX, double tCur, int targetDimension, boolean isSammon, double avgDist, double sumOfSquareDist)
        throws MPIException {

        final double [] sigmaValues = new double [ParallelOps.threadCount];
        IntStream.range(0, ParallelOps.threadCount).forEach(i -> sigmaValues[i] = 0.0);

        if (ParallelOps.threadCount > 1) {
            launchHabaneroApp(
                () -> forallChunked(
                    0, ParallelOps.threadCount - 1,
                    (threadIdx) -> sigmaValues[threadIdx] =
                        calculateStressInternal(threadIdx, preX, targetDimension, tCur, isSammon, avgDist)));
            // Sum across threads and accumulate to zeroth entry
            IntStream.range(1, ParallelOps.threadCount).forEach(
                i -> {
                    sigmaValues[0] += sigmaValues[i];
                });
        }
        else {
            sigmaValues[0] = calculateStressInternal(0, preX, targetDimension, tCur, isSammon, avgDist);
        }

        if (ParallelOps.procCount > 1) {
            sigmaValues[0] = ParallelOps.allReduce(sigmaValues[0]);
        }
        return sigmaValues[0] / sumOfSquareDist;
    }

    private static double calculateStressInternal(
        int threadIdx, double[][] preX, int targetDim, double tCur,
        boolean isSammon, double avgDist) {

        double sigma = 0.0;
        double diff = 0.0;
        if (tCur > 10E-10) {
            diff = Math.sqrt(2.0 * targetDim) * tCur;
        }

        int pointCount =
            ParallelOps.threadRowCounts[threadIdx] * ParallelOps.globalColCount;

        for (int i = 0; i < pointCount; ++i) {
            int procLocalPnum =
                i + ParallelOps.threadPointStartOffsets[threadIdx];
            double origD = distances.getValue(procLocalPnum);
            double weight = weights.getValue(procLocalPnum);

            if (origD < 0 || weight == 0) {
                continue;
            }

            weight = isSammon ? 1.0 / Math.max(origD, 0.001 * avgDist) : weight;
            int globalPointStart =
                procLocalPnum + ParallelOps.procPointStartOffset;
            int globalRow = globalPointStart / ParallelOps.globalColCount;
            int globalCol = globalPointStart % ParallelOps.globalColCount;

            double euclideanD = globalRow != globalCol ? calculateEuclideanDist(
                preX, targetDim, globalRow, globalCol) : 0.0;

            double heatD = origD - diff;
            double tmpD = origD >= diff ? heatD - euclideanD : 0;
            sigma += weight * tmpD * tmpD;
        }
        return sigma;
    }

    private static double calculateEuclideanDist(
        double[][] vectors, int targetDim, int i, int j) {
        double dist = 0.0;
        for (int k = 0; k < targetDim; k++) {
            double diff = vectors[i][k] - vectors[j][k];
            dist += diff * diff;
        }

        dist = Math.sqrt(dist);
        return dist;
    }

    static double[][] generateInitMapping(int numDataPoints,
                                          int targetDim) {
        double matX[][] = new double[numDataPoints][targetDim];
        // Use Random class for generating random initial mapping solution.
        Random rand = new Random(System.currentTimeMillis());
        for (int i = 0; i < numDataPoints; i++) {
            for (int j = 0; j < targetDim; j++) {
                if(rand.nextBoolean())
                    matX[i][j] = rand.nextDouble();
                else
                    matX[i][j] = -rand.nextDouble();
            }
        }
        return matX;
    }

    private static DoubleStatistics calculateStatistics(
        RefObj<Integer> missingDistCount)
        throws MPIException {
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
                            threadIdx, missingDistCounts)));
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
                0, missingDistCounts);
        }

        if (ParallelOps.procCount > 1) {
            threadDistanceSummaries[0] =
                ParallelOps.allReduce(threadDistanceSummaries[0]);
            missingDistCounts[0] = ParallelOps.allReduce(missingDistCounts[0]);
        }
        missingDistCount.setValue(missingDistCounts[0]);
        return threadDistanceSummaries[0];
    }

    private static void readDistancesAndWeights() {
        distances = BinaryReader.readRowRange(
            config.distanceMatrixFile, ParallelOps.procRowRange,
            ParallelOps.globalColCount, byteOrder, config.isMemoryMapped, true);
        if (config.distanceTransform != 1.0){
            distances = BinaryReader.transform(d -> d < 0 ? d : Math.pow(d, config.distanceTransform), distances);
        }

        weights = Strings.isNullOrEmpty(config.weightMatrixFile) ? BinaryReader
            .readConstant(1.0) : BinaryReader.readRowRange(
            config.weightMatrixFile, ParallelOps.procRowRange,
            ParallelOps.globalColCount, byteOrder, config.isMemoryMapped,
            false);

        /* TODO See if arrays would improve performance, especially within CG matrix multiply routine as it gets called many times*/
        distanceArray = new double[ParallelOps.procRowCount][ParallelOps.globalColCount];
        weightsArray = new double[ParallelOps.procRowCount][ParallelOps.globalColCount];

    }

    private static DoubleStatistics calculateStatisticsInternal(
        int threadIdx, int[] missingDistCounts) {

        DoubleStatistics stat = new DoubleStatistics();
        int pointCount =  ParallelOps.threadRowCounts[threadIdx] *
                          ParallelOps.globalColCount;
        for (int i = 0; i < pointCount; ++i){
            int procLocalPnum =
                i + ParallelOps.threadPointStartOffsets[threadIdx];
            double origD = distances.getValue(procLocalPnum);
            double weight = weights.getValue(procLocalPnum);
            if (origD < 0) {
                // Missing distance
                ++missingDistCounts[threadIdx];
                continue;
            }
            if (weight == 0) continue; // Ignore zero weights
            stat.accept(origD);
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
        byteOrder =
            config.isBigEndian ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN;
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
            System.out.println(e);
        }
        return Optional.fromNullable(null);
    }
}
