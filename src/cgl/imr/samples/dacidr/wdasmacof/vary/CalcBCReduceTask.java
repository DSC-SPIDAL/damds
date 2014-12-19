package cgl.imr.samples.dacidr.wdasmacof.vary;

import java.util.List;

import cgl.imr.base.Key;
import cgl.imr.base.ReduceOutputCollector;
import cgl.imr.base.ReduceTask;
import cgl.imr.base.TwisterException;
import cgl.imr.base.Value;
import cgl.imr.base.impl.JobConf;
import cgl.imr.base.impl.ReducerConf;
import cgl.imr.types.StringKey;

/**
 * @author Yang Ruan, yangruan@cs.indiana.edu
 * 
 *         This class simply combine the matrix blocks produce by the
 *         CalcBCMapTasks
 * 
 */
public class CalcBCReduceTask implements ReduceTask {

	int N = 0;

	public void configure(JobConf jobConf, ReducerConf mapConf)
			throws TwisterException {
		N = Integer.parseInt(jobConf.getProperty(DAMDS2.PROP_N));
	}

	/**
	 * Perform the actual operation of combining the set of matrix blocks into a
	 * single matrix. The offset comes with the data gives the location of the
	 * data block in the final matrix.
	 */
	public void reduce(ReduceOutputCollector collector, Key key,
			List<Value> values) throws TwisterException {

		// need to know the size of the centroids.
		// The reduce inputs are hard coded to be ByteValues.
		// TODO change this to generic values later.
		MDSMatrixData firstData = (MDSMatrixData) values.get(0);
		int width = firstData.getWidth();

		double results[][] = new double[N][width];
		double tmpData[][] = null;
		MDSMatrixData mData = null;
		int height = 0;
		int offset = 0;
		for (Value val : values) {
			mData = (MDSMatrixData) val;
			height = mData.getHeight();
			offset = mData.getRowOffset();
			tmpData = mData.getData();
			for (int m = 0; m < height; m++) {
				for (int n = 0; n < width; n++) {
					results[m + offset][n] = tmpData[m][n];
				}
			}
		}

		MDSMatrixData rowX = new MDSMatrixData(results, N, width);
		collector.collect(new StringKey("calc-x-reduce-combiner-key"), rowX);

	}

	public void close() throws TwisterException {
		// TODO Auto-generated method stub

	}

}
