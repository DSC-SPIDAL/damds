package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.damds.ParallelOps;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.concurrent.TimeUnit;

public class CGLoopTimings {
    public static enum TimingTask{
        MM, INNER_PROD_PAP, INNER_PROD_R
    }

    private Stopwatch timerMM = Stopwatch.createUnstarted();
    private Stopwatch timerInnerProdPAP = Stopwatch.createUnstarted();
    private Stopwatch timerInnerProdR = Stopwatch.createUnstarted();

    private long tMM;
    private long tInnerProdPAP;
    private long tInnerProdR;

    private long countMM;
    private long countInnerProdPAP;
    private long countInnerProdR;

    public void startTiming(TimingTask task){
        switch (task){
            case MM:
                timerMM.start();
                ++countMM;
                break;
            case INNER_PROD_PAP:
                timerInnerProdPAP.start();
                ++countInnerProdPAP;
                break;
            case INNER_PROD_R:
                timerInnerProdR.start();
                ++countInnerProdR;
                break;
        }
    }

    public void endTiming(TimingTask task){
        switch (task){
            case MM:
                timerMM.stop();
                tMM += timerMM.elapsed(TimeUnit.MILLISECONDS);
                timerMM.reset();
                break;
            case INNER_PROD_PAP:
                timerInnerProdPAP.stop();
                tInnerProdPAP += timerInnerProdPAP.elapsed(TimeUnit.MILLISECONDS);
                timerInnerProdPAP.reset();
                break;
            case INNER_PROD_R:
                timerInnerProdR.stop();
                tInnerProdR += timerInnerProdR.elapsed(TimeUnit.MILLISECONDS);
                timerInnerProdR.reset();
                break;
        }
    }

    public double getTotalTime(TimingTask task){
        switch (task){
            case MM:
                return tMM;
            case INNER_PROD_PAP:
                return tInnerProdPAP;
            case INNER_PROD_R:
                return tInnerProdR;
        }
        return  0.0;
    }

    public double getAverageTime(TimingTask task){
        switch (task){
            case MM:
                return tMM *1.0/ countMM;
            case INNER_PROD_PAP:
                return tInnerProdPAP *1.0/ countInnerProdPAP;
            case INNER_PROD_R:
                return tInnerProdR *1.0/ countInnerProdR;
        }
        return  0.0;
    }

    public long[] getCountDistribution(TimingTask task) throws
        MPIException {
        LongBuffer mpiOnlyTimingBuffer =  ParallelOps.mpiOnlyBuffer;
        mpiOnlyTimingBuffer.position(0);
        long [] mpiOnlyTimingArray = new long[ParallelOps.worldProcsCount];
        switch (task){
            case MM:
                mpiOnlyTimingBuffer.put(countMM);
                break;
            case INNER_PROD_PAP:
                mpiOnlyTimingBuffer.put(countInnerProdPAP);
                break;
            case INNER_PROD_R:
                mpiOnlyTimingBuffer.put(countInnerProdR);
                break;
        }
        ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
        return mpiOnlyTimingArray;
    }

}
