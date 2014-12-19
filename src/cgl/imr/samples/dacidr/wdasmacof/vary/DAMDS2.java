package cgl.imr.samples.dacidr.wdasmacof.vary;

/**
 * @author Yang Ruan: yangruan@cs.indiana.edu
 *		Twister-WDAMDS is implemented based on Twister-MDS.
 */

import cgl.imr.base.Key;
import cgl.imr.base.TwisterException;
import cgl.imr.base.TwisterModel;
import cgl.imr.base.TwisterMonitor;
import cgl.imr.base.impl.GenericCombiner;
import cgl.imr.base.impl.JobConf;
import cgl.imr.client.TwisterDriver;
import cgl.imr.types.DoubleArray;
import cgl.imr.types.DoubleValue;
import cgl.imr.types.StringValue;
import org.safehaus.uuid.UUIDGenerator;

import java.io.*;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public class DAMDS2 {

	private static UUIDGenerator uuidGen = UUIDGenerator.getInstance();
	public static final String DATE_FORMAT_NOW = "yyyy-MM-dd";
	public static String PROP_BZ = "mat_mult_block_size";
	public static String PROP_N = "prop_N";
	public static String PROP_TCUR = "prop_T_cur";
	public static String PROP_D = "prop_target_dim";
	public static String PROP_ALPHA = "prop_alpha";
	public static String PROP_SAMMON = "prop_sammon";
	public static String PROP_AVG_D = "prop_avg_d";
	private static double SEQUENTIAL_TIME = 0;
	
	public static int BLOCK_SIZE = 64;

	public static String date() {
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
		return sdf.format(cal.getTime());
	}

	/*
	 * This is the maximum determinant value for using Inverse or pseudo
	 * inverse.
	 */
	public static double MAX_DET = 10000000;
	public static int MAX_ITER = 10000;

	static int[] ids;
	static int N;
	static int D;
	static double avgOrigDistance;
	static double sumOrigDistanceSquare;
	static double avgOrigDistanceSquare;
	static double maxOrigDistance;
	static double tMax;
	static double tCur = 0.0;
	static double tMin;
	static double alpha;
	static int CG_ITER;
	static double CG_REAL_ITER = 0;
	static double SMACOF_REAL_ITER = 0;
	static double CG_THRESHOLD = 1;
    static boolean sammonMapping = false;

	public static void main(String[] args) {

		double beginTime = System.currentTimeMillis();

		if (args.length != 15) {
			System.out.println("Usage: ");
			System.out.println("[1. Num map tasks ]");
			System.out.println("[2. Input Folder]");
			System.out.println("[3. Input File Prefix]");
			System.out.println("[4. Input Weight Prefix]");
			System.out.println("[5. IDs File ]");
			System.out.println("[6. Label Data File ]");
			System.out.println("[7. Output File ]");
			System.out.println("[8. Threshold value ]");
			System.out.println("[9. The Target Dimension ]");
			System.out.println("[10. Cooling parameter (alpha) ]");
			System.out.println("[11. Input Data Size]");
			System.out.println("[12. Final Weight Prefix]");
			System.out.println("[13. CG iteration num]");
			System.out.println("[14. CG Error Threshold]");
			System.out.println("[15. Sammon mapping (boolean)]");
			System.exit(0);
		}

		int numMapTasks = Integer.valueOf(args[0]);
		String inputFolder = args[1];
		String inputPrefix = args[2];
		String weightPrefix = args[3];
		String idsFile = args[4];
		String labelsFile = args[5];
		String outputFile = args[6];
		double threshold = Double.valueOf(args[7]);
		D = Integer.valueOf(args[8]);
		alpha = Double.valueOf(args[9]);
		N = Integer.parseInt(args[10]);
		String finalWeightPrefix = args[11];
		CG_ITER = Integer.parseInt(args[12]);
		CG_THRESHOLD = Double.parseDouble(args[13]);
        sammonMapping = Boolean.parseBoolean(args[14]);

		System.out.println("[1. Num map tasks ]:\t" + numMapTasks);
		System.out.println("[2. Input Folder]:\t" + inputFolder);
		System.out.println("[3. Input File Prefix]:\t" + inputPrefix);
		System.out.println("[4. Weighted File Prefix]:\t" + weightPrefix);
		System.out.println("[5. IDs File ]:\t" + idsFile);
		System.out.println("[6. Label Data File ]:\t" + labelsFile);
		System.out.println("[7. Output File ]:\t" + outputFile);
		System.out.println("[8. Threshold value ]:\t" + threshold);
		System.out.println("[9. The Target Dimension ]:\t" + D);
		System.out.println("[10. Cooling parameter (alpha) ]:\t" + alpha);
		System.out.println("[11. Input Data Size]:\t" + N);
		System.out.println("[12. Final Weight Prefix]:\t" + finalWeightPrefix);
		System.out.println("[13. CG Iterations]:\t" + CG_ITER);
		System.out.println("[14. CG Threshold]:\t" + CG_THRESHOLD);
        System.out.println("[15. Sammon mapping]:\t" +sammonMapping);
		
		try {
			performInitialCalculations(numMapTasks, inputFolder, inputPrefix, 
					idsFile, weightPrefix);
			System.out.println(" N: " + N);
			System.out.println(" AvgOrgDistance: " + avgOrigDistance);
			System.out.println(" SumSquareOrgDistance: " + sumOrigDistanceSquare);
			System.out.println(" MaxOrigDistance: " + maxOrigDistance);

			double[][] preX = generateInitMapping(N, D);

			TwisterModel stressDriver = configureCalculateStress(numMapTasks,
					inputFolder, inputPrefix, weightPrefix, idsFile);
			Double stress = null;

			Double preStress = calculateStress(stressDriver, preX, numMapTasks);
			System.out.println("Initial Stress: " + preStress);

			double diffStress = 10 * threshold;  //starting value
			int iter = 0;

			double X[][] = null;
			double BC[][] = null;

			// Configuring BC MapReduce driver.
			TwisterModel bcDriver = 
					configureCalculateBC(numMapTasks, inputFolder, 
							inputPrefix, weightPrefix, idsFile);
			
			TwisterModel mmDriver = 
					configureMatrixMutiply(numMapTasks, inputFolder, inputPrefix, weightPrefix, idsFile);

			double QoR1 = 0;
			double QoR2 = 0;

			double avgOrigDist = avgOrigDistance;

			
			tMax = calculateMaxT(maxOrigDistance, D);
			tMin = (0.01 * tMax < 0.01) ? 0.01 * tMax : 0.01;
			tCur = alpha * tMax;
			
			double endTime = System.currentTimeMillis();
			System.out.println("Upto the loop took =" + (endTime - beginTime)
					/ 1000 + " Seconds.");

			iter = 0;
			
			while (true) {
				preStress = calculateStress(stressDriver, preX, numMapTasks);
				diffStress = threshold + 1.0;
				
				System.out.println("###############################");
				System.out.printf("# T_Cur = %.10g\n", tCur);
				System.out.println("###############################");
				
				while ( diffStress >= threshold ) {
					BC = calculateBC(bcDriver, preX);
					
					X = conjugateGradient(mmDriver, BC, preX);
					
					stress = calculateStress(stressDriver, X, numMapTasks);
					diffStress = preStress - stress;
					preStress = stress;
					preX = MatrixUtils.copy(X);

					iter++;
					
					if ((iter % 10 == 0) || (iter >= MAX_ITER)) {
						System.out.println("Iteration ## " + iter + " completed. " + threshold + " " + diffStress + " " + stress);
					}
					++SMACOF_REAL_ITER;
				}
				
				System.out.println("Iteration ## " + iter + " completed. " + threshold + " " + diffStress + " " + stress);
				System.out.println();
				
				if (tCur == 0)
					break;
				tCur *= alpha;
				if (tCur < tMin)
					tCur = 0;
				iter = 0;
			}
			
			QoR1 = stress / (N * (N - 1) / 2);
			QoR2 = QoR1 / (avgOrigDist * avgOrigDist);

			System.out.println("Normalize1 = " + QoR1 + " Normalize2 = " + QoR2);
			System.out.println("Average of Delta(original distance) = "	+ avgOrigDist);

			endTime = System.currentTimeMillis();

			if (labelsFile.endsWith("NoLabel")) {
				writeOuput(X, outputFile);
			} else {
				writeOuput(X, labelsFile, outputFile);
			}
			bcDriver.close();			
			mmDriver.close();
			stressDriver.close();
			
			
			TwisterModel finalStressDriver = configureCalculateStress(numMapTasks,
					inputFolder, inputPrefix, finalWeightPrefix, idsFile);

			Double finalStress = calculateStress(finalStressDriver, X, numMapTasks);

			System.out
					.println("===================================================");
			System.out.println("CG REAL ITER:" + CG_REAL_ITER);
			System.out.println("SMACOF REAL ITER: " + SMACOF_REAL_ITER);
			System.out.println("For CG iter: " + CG_REAL_ITER / SMACOF_REAL_ITER + "\tFinal Result is:\t" + finalStress + "\t" + (endTime - beginTime) / 1000);
			System.out.println("Sequential Computation Time: " + SEQUENTIAL_TIME);
			System.out
					.println("===================================================");

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		System.exit(0);
	}
	
	private static double[][] conjugateGradient(TwisterModel mmDriver, 
			double[][] BC, double[][] preX) throws TwisterException{
		double[][] X = null;
		double[][] r = new double[N][D];
		double[][] p = new double[N][D];
		
		X = preX;
		r = calcMM(mmDriver, X);

		for(int i = 0; i < N; ++i)
			for(int j = 0; j < D; ++j){
				p[i][j] = BC[i][j] - r[i][j];
				r[i][j] = p[i][j];
			}
		
		int cgCount = 0;
		double rTr = innerProductCalculation(r);

		//System.out.println("1");
		while(cgCount < CG_ITER){
			cgCount++;
			++CG_REAL_ITER;
			//System.out.println("2");
			//calculate alpha
			double[][] Ap = calcMM(mmDriver, p);
			
			double alpha = rTr
					/innerProductCalculation(p, Ap);

			//update Xi to Xi+1
			for(int i = 0; i < N; ++i)
				for(int j = 0; j < D; ++j)
					X[i][j] += alpha * p[i][j];

			if (rTr < CG_THRESHOLD) {
				break;
			}
			
			//update ri to ri+1
			for(int i = 0; i < N; ++i)
				for(int j = 0; j < D; ++j)
					r[i][j] = r[i][j] - alpha * Ap[i][j];
			
			//calculate beta
			double rTr1 = innerProductCalculation(r);
			double beta = rTr1/rTr;
			rTr = rTr1;
			
			//update pi to pi+1
			for(int i = 0; i < N; ++i)
				for(int j = 0; j < D; ++j)
					p[i][j] = r[i][j] + beta * p[i][j];

			
		}
		return X;
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
				System.out.println("ERROR: Invalid lable");
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
			writer.println(labels.get(String.valueOf(ids[i]))); // print label
																// value, which
																// is
			// ONE for all data.
		}
		writer.flush();
		writer.close();
	}

	private static double[][] calculateBC(TwisterModel bcDriver, double[][] preX)
			throws TwisterException {

		MDSMatrixData BCMatData = null;
		MDSMatrixData preXMatData = new MDSMatrixData(preX, N, preX[0].length);
		preXMatData.setCurT(tCur);
		
		String memCacheKey = bcDriver.addToMemCache(preXMatData);
		TwisterMonitor monitor = bcDriver.runMapReduceBCast(new StringValue(
				memCacheKey));
		monitor.monitorTillCompletion();
		SEQUENTIAL_TIME += monitor.getTotalSequentialTimeSeconds();
		bcDriver.cleanMemCache(memCacheKey);
		GenericCombiner combiner = (GenericCombiner) bcDriver
				.getCurrentCombiner();
		if (!combiner.getResults().isEmpty()) {
			Key key = combiner.getResults().keySet().iterator().next();
			BCMatData = (MDSMatrixData) (combiner.getResults().get(key));
			return BCMatData.getData();
		} else {			
		/*
			if(status.isSuccess()){
				System.err.println("Combiner did not return any values. Terminating Job");
				bcDriver.close();
				System.exit(-1);
			}
			*/
			return null;
		}	

	}

	private static TwisterModel configureCalculateBC(int numMapTasks,
			String inputFolder, String inputPrefix, String weightPrefix, String idsFile) throws TwisterException {
		String jobID = "calc-BC-" + uuidGen.generateRandomBasedUUID();
		// we need only one reduce task to aggregate the parts of BC
		int numReducers = 1;

		// JobConfigurations
		JobConf jobConf = new JobConf(jobID);
		jobConf.setMapperClass(CalcBCMapTask.class);
		jobConf.setReducerClass(CalcBCReduceTask.class);
		jobConf.setCombinerClass(GenericCombiner.class);
		jobConf.setNumMapTasks(numMapTasks);
		jobConf.setNumReduceTasks(numReducers);
		jobConf.addProperty(PROP_BZ, String.valueOf(BLOCK_SIZE));
		jobConf.addProperty(PROP_N, String.valueOf(N));
		
		jobConf.addProperty(PROP_TCUR, String.valueOf(tCur));
		jobConf.addProperty(PROP_D, String.valueOf(D));
        jobConf.addProperty(PROP_SAMMON, String.valueOf(sammonMapping));
        jobConf.addProperty(PROP_AVG_D, String.valueOf(avgOrigDistance));
		jobConf.addProperty("InputFolder", inputFolder);
		jobConf.addProperty("InputPrefix", inputPrefix);
		jobConf.addProperty("WeightPrefix", weightPrefix);
		jobConf.addProperty("IdsFile", idsFile);

		TwisterModel bcDriver = new TwisterDriver(jobConf);
		bcDriver.configureMaps();
		return bcDriver;

	}

	private static double[][] calcMM(TwisterModel xDriver, double[][] X)
			throws TwisterException {

		MDSMatrixData rMatData = null;
		MDSMatrixData XMatData = new MDSMatrixData(X, X.length, X[0].length);
		String memCacheKey = xDriver.addToMemCache(XMatData);
		TwisterMonitor monitor = xDriver.runMapReduceBCast(new StringValue(
				memCacheKey));
		monitor.monitorTillCompletion();
		SEQUENTIAL_TIME += monitor.getTotalSequentialTimeSeconds();
		xDriver.cleanMemCache(memCacheKey);

		GenericCombiner combiner = (GenericCombiner) xDriver
				.getCurrentCombiner();
		if (!combiner.getResults().isEmpty()) {
			Key key = combiner.getResults().keySet().iterator().next();
			rMatData = (MDSMatrixData) (combiner.getResults().get(key));
			return rMatData.getData();
		} else {
			/*
			if(status.isSuccess()){
				System.err.println("Combiner did not return any values. Terminating Job");
				xDriver.close();
				System.exit(-1);
			}
			*/
			return null;
		}
		

	}

	private static TwisterModel configureMatrixMutiply(int numMapTasks,
			String inputFolder, String inputPrefix, String weightPrefix, String idsFile) throws TwisterException {
		String jobID = "calc-CG-" + uuidGen.generateRandomBasedUUID();
		// we need only one reduce task to aggregate the parts of X.
		int numReducers = 1;

		// JobConfigurations
		JobConf jobConf = new JobConf(jobID);
		jobConf.setMapperClass(MatrixMultiplyMapTask.class);
		jobConf.setReducerClass(MatrixMultiplyReduceTask.class);
		jobConf.setCombinerClass(GenericCombiner.class);
		jobConf.setNumMapTasks(numMapTasks);
		jobConf.setNumReduceTasks(numReducers);
		jobConf.addProperty(PROP_BZ, String.valueOf(BLOCK_SIZE));
		jobConf.addProperty(PROP_N, String.valueOf(N));
        jobConf.addProperty(PROP_SAMMON, String.valueOf(sammonMapping));
        jobConf.addProperty(PROP_AVG_D, String.valueOf(avgOrigDistance));
		jobConf.addProperty("InputFolder", inputFolder);
		jobConf.addProperty("InputPrefix", inputPrefix);
		jobConf.addProperty("WeightPrefix", weightPrefix);
		jobConf.addProperty("IdsFile", idsFile);

		TwisterModel xDriver = new TwisterDriver(jobConf);
		xDriver.configureMaps();
		return xDriver;
	}

	public static Double calculateStress(TwisterModel stressDriver,
			double[][] preX, int numMapTasks) throws TwisterException {
		double stress = 0;

		MDSMatrixData preXMatData = new MDSMatrixData(preX, N, preX[0].length);
		preXMatData.setCurT(tCur);
		
		String memCacheKey = stressDriver.addToMemCache(preXMatData);
		TwisterMonitor monitor = stressDriver
				.runMapReduceBCast(new StringValue(memCacheKey));
		monitor.monitorTillCompletion();
		SEQUENTIAL_TIME += monitor.getTotalSequentialTimeSeconds();
		stressDriver.cleanMemCache(memCacheKey);

		GenericCombiner combiner = (GenericCombiner) stressDriver
				.getCurrentCombiner();
		if (!combiner.getResults().isEmpty()) {
			Key key = combiner.getResults().keySet().iterator().next();
			stress = ((DoubleValue) combiner.getResults().get(key)).getVal();
			return new Double(stress / (sumOrigDistanceSquare));
		} else {
			/*
			if(status.isSuccess()){
				System.err.println("Combiner did not return any values. Terminating Job");
				stressDriver.close();
				System.exit(-1);
			}
			*/
			return null;
		}
		/*
		 * Stress divided by half of the sum of the square of the original
		 * distance to normalize it. We need it to be by half of the sum because
		 * when we calculate the stress we do it only for the upper triangular
		 * matrix. Otherwise we could just divide by the sum.
		 */
		// return stress/(sumOrigDistanceSquare/2);
		
	}

	public static TwisterModel configureCalculateStress(int numMapTasks,
			String inputFolder, String inputPrefix, String weightPrefix, String idsFile) throws TwisterException {
		String jobID = "stress-calc-" + uuidGen.generateRandomBasedUUID();
		// we need only one reducer for the above algorithm.
		int numReducers = 1;

		// JobConfigurations
		JobConf jobConf = new JobConf(jobID);
		jobConf.setMapperClass(StressMapTask.class);
		jobConf.setReducerClass(StressReduceTask.class);
		jobConf.setCombinerClass(GenericCombiner.class);
		jobConf.setNumMapTasks(numMapTasks);
		jobConf.setNumReduceTasks(numReducers);

		jobConf.addProperty(PROP_TCUR, String.valueOf(tCur));
		jobConf.addProperty(PROP_D, String.valueOf(D));
        jobConf.addProperty(PROP_SAMMON, String.valueOf(sammonMapping));
        jobConf.addProperty(PROP_AVG_D, String.valueOf(avgOrigDistance));
		jobConf.addProperty("InputFolder", inputFolder);
		jobConf.addProperty("InputPrefix", inputPrefix);
		jobConf.addProperty("WeightPrefix", weightPrefix);
		jobConf.addProperty("IdsFile", idsFile);
//		jobConf.setFaultTolerance();

		TwisterModel stressDriver = new TwisterDriver(jobConf);
		stressDriver.configureMaps();

		return stressDriver;
	}

	/**
	 * This method will be used to generate initial mapped data X(0), when there
	 * is no initial mapped data for the problem. Normally, target space
	 * 2-dimension space.
	 * 
	 * @param numDataPoints
	 * @param targetDim
	 * @return
	 */
	static double[][] generateInitMapping(int numDataPoints,
			int targetDim) {
		double matX[][] = new double[numDataPoints][targetDim];
		// Use Random class for generating random initial mapping solution.
		// For the test, set the Random seed in order to produce the same result
		// for the same problem.
		// Random rand = new Random(47);
		Random rand = new Random(System.currentTimeMillis()); // Real random
		// seed.
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


	private static double calculateMaxT(double maxOrigDistance, int targetDim) {
		double divider = Math.sqrt(2.0 * targetDim);

		return maxOrigDistance / divider;
	}
	
	public static double innerProductCalculation(double[][] a){
		int row = a.length, col = a[0].length;
		double sum = 0;
		for(int i = 0; i < row; ++i)
			for(int j = 0; j < col; ++j)
				sum += a[i][j] * a[i][j];
		return sum;
	}
	
	public static double innerProductCalculation(double[][] a, double[][] b){
		int row = a.length, col = a[0].length;
		double sum = 0;
		for(int i = 0; i < row; ++i)
			for(int j = 0; j < col; ++j)
				sum += a[i][j] * b[i][j];
		return sum;
	}
	
	/**
	 * This will set values for the following parameters. (i) N - size of the
	 * matrix, (ii)ids -ID vector, and (iii) averageOriginalDisance - average of
	 * sum of all the distances.
	 * 
	 * @param numMapTasks
	 *            - Number of map tasks.
	 * @param idsFile
	 *            - File containing IDs.
	 * @throws Exception
	 */
	private static void performInitialCalculations(int numMapTasks,
			String inputFolder, String inputPrefix, String idsFile, 
			String weightedPrefix) throws Exception {
			
		// JobConfigurations for calculating original average distance.
		JobConf jobConf = new JobConf("avg-original-distance"
				+ uuidGen.generateTimeBasedUUID());
		jobConf.setMapperClass(AvgOrigDistanceMapTask.class);
		jobConf.setReducerClass(AvgOrigDistanceReduceTask.class);
		jobConf.setCombinerClass(GenericCombiner.class);
		jobConf.setNumMapTasks(numMapTasks);
		jobConf.setNumReduceTasks(1);// One reducer is enough since we just
		// summing some numbers.
        jobConf.addProperty(PROP_SAMMON, String.valueOf(sammonMapping));
		jobConf.addProperty("InputFolder", inputFolder);
		jobConf.addProperty("InputPrefix", inputPrefix);
		jobConf.addProperty("WeightPrefix", weightedPrefix);
		jobConf.addProperty("IdsFile", idsFile);

		TwisterModel driver = null;
		TwisterMonitor monitor = null;
		GenericCombiner combiner;
		try {
			driver = new TwisterDriver(jobConf);
			driver.configureMaps();

			monitor = driver.runMapReduce();
			monitor.monitorTillCompletion();
			combiner = (GenericCombiner) driver.getCurrentCombiner();
			if (!combiner.getResults().isEmpty()) {
				Key key = combiner.getResults().keySet().iterator().next();
				avgOrigDistance = ((DoubleArray) combiner.getResults().get(key))
						.getData()[0];
				sumOrigDistanceSquare = ((DoubleArray) combiner.getResults()
						.get(key)).getData()[1];
				maxOrigDistance = ((DoubleArray) combiner.getResults().get(key)).getData()[2];

				// Ignoring the diagonal zeros from the average.
				double div = Math.pow(N, 2)-N;
				avgOrigDistance = avgOrigDistance/div;
				avgOrigDistanceSquare = sumOrigDistanceSquare / div;
			} else {
				System.err.println("Combiner did not return any values.");
				driver.close();
				System.exit(-1);
			}

		} catch (TwisterException e) {
			driver.close();
			throw e;
		}
		driver.close();
	}
}
