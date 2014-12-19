package cgl.imr.samples.dacidr.wdasmacof.vary;

import java.util.List;

import cgl.imr.base.Key;
import cgl.imr.base.ReduceOutputCollector;
import cgl.imr.base.ReduceTask;
import cgl.imr.base.TwisterException;
import cgl.imr.base.Value;
import cgl.imr.base.impl.JobConf;
import cgl.imr.base.impl.ReducerConf;
import cgl.imr.types.DoubleValue;
import cgl.imr.types.StringKey;

/**
 * @author Yang Ruan, yangruan@cs.indiana.edu
 * 
 *         Simply add the partial stress and make it a total. Lot of stesses :)
 * 
 */
public class StressReduceTask implements ReduceTask {

	JobConf jobConf;

	public void reduce(ReduceOutputCollector collector, Key key,
			List<Value> values) throws TwisterException {

		double totalStress = 0;
		for (Value val : values) {
			DoubleValue sigmaVal = (DoubleValue) val;
			totalStress += sigmaVal.getVal();

		}
		// Send the total stress as string
		collector.collect(new StringKey("stress-map-to-combine-key"),
				new DoubleValue(totalStress));

	}

	@Override
	public void configure(JobConf jobConf, ReducerConf reducerConf)
			throws TwisterException {
		// TODO Auto-generated method stub

	}

	@Override
	public void close() throws TwisterException {
		// TODO Auto-generated method stub

	}

}
