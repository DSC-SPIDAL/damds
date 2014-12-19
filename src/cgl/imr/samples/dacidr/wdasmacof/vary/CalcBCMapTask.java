package cgl.imr.samples.dacidr.wdasmacof.vary;

import cgl.imr.base.*;
import cgl.imr.base.impl.JobConf;
import cgl.imr.base.impl.MapperConf;
import cgl.imr.types.StringKey;
import cgl.imr.types.StringValue;
import cgl.imr.worker.MemCache;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * 
 * @author Yang Ruan, yangruan@cs.indiana.edu Some code segments contain
 *         in this class are directly inherited from the C# version of the
 *         shared memory MDS program wirtten by my collegue Seung-Hee Bae and Jaliya
 * 
 *         This class performs the partial calculation of the following matrrix
 *         multiplication. BC= BofZ * preX. Please see ParallelMDSMapReduce.java
 *         for a complete description of the terms and the algorithm.
 * 
 * @author Seung-Hee Bae: sebae@cs.indiana.edu
 *			Modify Twister-MDS version of this class for Twister-DAMDS version.
 */

public class CalcBCMapTask implements MapTask {
    private boolean sammonMapping = false;
    private double averageOriginalDistance = 0.0;
	int mapNo = 0;
	private int bz = 0; // this is to hold the block size of the block matrix
	private int N = 0;

	private double tCur = 0.0;
	private int targetDim = 3;
	
	MDSShortMatrixData deltaMatData;
	int blockOffset = 0;
	int blockHeight = 0;

	private short[][] deltaBlock = null;
	private short[][] weights = null;

	float[][] BofZ = null;
	private JobConf jobConf;

	public void map(MapOutputCollector collector, Key key, Value val)
			throws TwisterException {

		StringValue memCacheKey = (StringValue) val;
		MDSMatrixData mData = (MDSMatrixData) (MemCache.getInstance().get(
				jobConf.getJobId(), memCacheKey.toString()));

		double[][] preX = mData.getData();
		double tmpCurT = mData.getCurT();
		
		if (tmpCurT != tCur) {
			tCur = tmpCurT;
		}

		// Now calculate the BofZ block using the bocks of matrices we have.
		// Each map task calculates the m th block of the matrix, where m is the
		// map task number.
		calculateBofZ(preX);

		// Next we can calculate the BofZ * preX.
		double[][] C = MatrixUtils.matrixMultiply(BofZ, preX, blockHeight,
				preX[0].length, N, bz);

		// Send C with the map task number to a reduce task. Which will simply
		// combine
		// these parts and form the N x d matrix.
		MDSMatrixData newMData = new MDSMatrixData(C, blockHeight, C[0].length,
				mData.getRow(), blockOffset);
		collector.collect(new StringKey("BC-calc-to-reduce-key"), newMData);

	}

	/**
	 * Calculation of partial BofZ matrix block.
	 * 
	 */
	private void calculateBofZ(double[][] preX) {
		int tmpI = 0;
		
		double vBlockValue = (double) -1;

		double diff = 0;
		if (tCur > 10E-10) {
			diff = Math.sqrt(2.0 * targetDim)  * tCur;
		}
        double weightMultiply = sammonMapping ? 1.0/Short.MAX_VALUE : 1.0;
		for (int i = blockOffset; i < blockOffset + blockHeight; i++) {
			tmpI = i - blockOffset;
			BofZ[tmpI][i] = 0;
			for (int j = 0; j < N; j++) {
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
				// separately.
				if (i != j) {
                    double weight = weights[tmpI][j] * weightMultiply;
                    if(weight != 0){
						double dist = calculateDistance(preX, preX[0].length, i, j);
						double origD = deltaBlock[tmpI][j] / (double)Short.MAX_VALUE;
						if (dist >= 1.0E-10 && diff < origD) {
							BofZ[tmpI][j] = (float) (weight * vBlockValue * (origD - diff) / dist);
						} else {
							BofZ[tmpI][j] = 0;
						}
						
						BofZ[tmpI][i] += -BofZ[tmpI][j];
					}
				}
			}
		}
	}

	/**
	 * Simple distance calculation
	 * 
	 * @param origMat
	 *            - original matrix
	 * @param vecLength
	 *            - width of the matrix.
	 * @param i
	 *            , j means the positions for the distance calculation.
	 * @return
	 */
	public static double calculateDistance(double[][] origMat, int vecLength,
			int i, int j) {
		/*
		 * i and j is the index of first dimension, actually the index of each
		 * points. the length of the second dimension is the length of the
		 * vector.
		 */
		double dist = 0;
		for (int k = 0; k < vecLength; k++) {
			double diff = origMat[i][k] - origMat[j][k];
			dist += diff * diff;
		}

		dist = Math.sqrt(dist);
		return dist;
	}

	/**
	 * During this configuration step, the map task will load two matrix blocks.
	 * 1. Block from the delta matrix 2. Block from V matrix. The data files
	 * contains the matrix block data in binary form.
	 */
	public void configure(JobConf jobConf, MapperConf mapConf)
			throws TwisterException {
        sammonMapping = Boolean.parseBoolean(jobConf.getProperty(DAMDS2.PROP_SAMMON));
        averageOriginalDistance = Double.parseDouble(jobConf.getProperty(DAMDS2.PROP_AVG_D));
		String inputFolder = jobConf.getProperty("InputFolder");
		String inputPrefix = jobConf.getProperty("InputPrefix");
		String weightPrefix = jobConf.getProperty("WeightPrefix");
		String fileName = (inputFolder + "/" + inputPrefix + mapConf.getMapTaskNo())
				.replaceAll("//", "/");
		String weightName = (inputFolder + "/" + weightPrefix + mapConf.getMapTaskNo())
				.replaceAll("//", "/");
		
		String idsFile = jobConf.getProperty("IdsFile");

		try {
			BufferedReader br = new BufferedReader(new FileReader(idsFile));
			String line;
			String[] tokens;
			deltaMatData = new MDSShortMatrixData();
			while((line = br.readLine())!=null){
				tokens = line.split("\t");
				if(Integer.parseInt(tokens[0]) == mapConf.getMapTaskNo()){
					deltaMatData.setHeight(Integer.parseInt(tokens[1]));
					deltaMatData.setWidth(Integer.parseInt(tokens[2]));
					deltaMatData.setRow(Integer.parseInt(tokens[3]));
					deltaMatData.setRowOFfset(Integer.parseInt(tokens[4]));
				}
			}
			br.close();
		}
		catch(IOException e){
			e.printStackTrace();
		}
		//System.out.println(mapConf.getMapTaskNo() + " " + fileName);
		try {
			deltaMatData.loadDeltaFromBinFile(fileName);
            weights = sammonMapping ? FileOperation.loadSammonWeights(deltaMatData.data, averageOriginalDistance, deltaMatData.getHeight(),
                    deltaMatData.getWidth()) : FileOperation.loadWeights(weightName, deltaMatData.getHeight(),
                    deltaMatData.getWidth());
		} catch (Exception e) {
			throw new TwisterException(e);
		}			

		this.jobConf = jobConf;

		bz = Integer.parseInt(jobConf.getProperty(DAMDS2.PROP_BZ));
		tCur = Double.parseDouble(jobConf.getProperty(DAMDS2.PROP_TCUR));
		targetDim = Integer.parseInt(jobConf.getProperty(DAMDS2.PROP_D));
		
		deltaBlock = deltaMatData.getData();
		
		blockOffset = deltaMatData.getRowOffset();
		blockHeight = deltaMatData.getHeight();
		N = deltaMatData.getWidth();
		BofZ = new float[blockHeight][N];
		mapNo = mapConf.getMapTaskNo();
	}

	public void close() throws TwisterException {
		// TODO Auto-generated method stub
	}
}
