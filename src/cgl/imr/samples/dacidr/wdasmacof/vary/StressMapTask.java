package cgl.imr.samples.dacidr.wdasmacof.vary;

import cgl.imr.base.*;
import cgl.imr.base.impl.JobConf;
import cgl.imr.base.impl.MapperConf;
import cgl.imr.types.DoubleValue;
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
 *         shared memory MDS program wirtten by my collegue Seung-Hee Bae.
 * 
 *         This class performs the partial calculation of the stress between two
 *         vectors.
 *
 * @author Seung-Hee Bae: sebae@cs.indiana.edu
 *			Modify Twister-MDS version for implementing Twister-DAMDS by Seung-Hee.
 */
public class StressMapTask implements MapTask {

    private boolean sammonMapping = false;
    private double averageOriginalDistance = 0.0;
	private double tCur = 0.0;
	private int targetDim = 3;
	private int rowHeight;
	private int N;

	MDSShortMatrixData deltaBlock;
	short[][] weights;
	MapperConf mapConf;
	JobConf jobConf;

	@Override
	public void configure(JobConf jobConf, MapperConf mapConf)
			throws TwisterException {
		this.mapConf = mapConf;
		this.jobConf = jobConf;

		tCur = Double.parseDouble(jobConf.getProperty(DAMDS2.PROP_TCUR));
		targetDim = Integer.parseInt(jobConf.getProperty(DAMDS2.PROP_D));

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
		try{
			BufferedReader br = new BufferedReader(new FileReader(idsFile));
			String line;
			String[] tokens;
			deltaBlock = new MDSShortMatrixData();
			while((line = br.readLine())!=null){
				tokens = line.split("\t");
				if(Integer.parseInt(tokens[0]) == mapConf.getMapTaskNo()){
					deltaBlock.setHeight(Integer.parseInt(tokens[1]));
					deltaBlock.setWidth(Integer.parseInt(tokens[2]));
					deltaBlock.setRow(Integer.parseInt(tokens[3]));
					deltaBlock.setRowOFfset(Integer.parseInt(tokens[4]));
				}
			}
			br.close();
		}
		catch(IOException e) {
			e.printStackTrace();
		}

		try {
			deltaBlock.loadDeltaFromBinFile(fileName);
            weights = sammonMapping ? FileOperation.loadSammonWeights(deltaBlock.data, averageOriginalDistance,
                    deltaBlock.getHeight(),
                    deltaBlock.getWidth()) : FileOperation.loadWeights(weightName, deltaBlock.getHeight(),
                    deltaBlock.getWidth());
	
			
			rowHeight = deltaBlock.getHeight();
			N = deltaBlock.getWidth();
			
		} catch (Exception e) {
			throw new TwisterException(e);
		}
	}

	public void map(MapOutputCollector collector, Key key, Value val)
			throws TwisterException {

		StringValue memCacheKey = (StringValue) val;
		MDSMatrixData mData = (MDSMatrixData) (MemCache.getInstance().get(
				jobConf.getJobId(), memCacheKey.toString()));

		int rowOffset = deltaBlock.getRowOffset();

		double[][] preXData = mData.getData();
		double tmpCurT = mData.getCurT();
		
		short deltaMatData[][] = deltaBlock.getData();

		int tmpI = 0;

		double sigma = 0.0;
		
		if (tmpCurT != tCur) {
			tCur = tmpCurT;
		}
		
		double diff = 0;
		if (tCur > 10E-10) {
			diff = Math.sqrt(2.0 * targetDim) * tCur;
		}
        double weightMultiply = sammonMapping ? 1.0/Short.MAX_VALUE : 1.0;
		for (int i = rowOffset; i < rowOffset + rowHeight; i++) {
			tmpI = i - rowOffset;
			for (int j = 0; j < N; j++) {
                double weight = weights[tmpI][j] * weightMultiply;
                if(weight != 0){
					double dist;
					if (j != i) {
						dist = calculateDistance(preXData, preXData[0].length, i, j);
					} else {
						dist = 0;
					}
					double origD = deltaMatData[tmpI][j] / (double)Short.MAX_VALUE;
					double heatDist = origD - diff;
					double d = origD >= diff 
								? heatDist - dist : 0;
					sigma += weight * d * d;
				}
			}
		}
		// Send the partial sigma.
		collector.collect(new StringKey("stress-map-to-reduce-key"),
				new DoubleValue(sigma));
	}

	private static double calculateDistance(double[][] origMat, int vecLength,
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

	public void close() throws TwisterException {
		// TODO Auto-generated method stub
	}
}