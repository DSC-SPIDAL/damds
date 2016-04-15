package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.damds.ParallelOps;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.concurrent.TimeUnit;

public class StressLoopTimings {
    public static enum TimingTask{
        BC,CG, STRESS
    }

    private Stopwatch timerBC = Stopwatch.createUnstarted();
    private Stopwatch timerCG = Stopwatch.createUnstarted();
    private Stopwatch timerStress = Stopwatch.createUnstarted();

    private long tBC;
    private long tCG;
    private long tStress;

    private long countBC;
    private long countCG;
    private long countStress;

    public void startTiming(TimingTask task){
        switch (task){
            case BC:
                timerBC.start();
                ++countBC;
                break;
            case CG:
                timerCG.start();
                ++countCG;
                break;
            case STRESS:
                timerStress.start();
                ++countStress;
                break;
        }
    }

    public void endTiming(TimingTask task){
        switch (task){
            case BC:
                timerBC.stop();
                tBC += timerBC.elapsed(TimeUnit.MILLISECONDS);
                timerBC.reset();
                break;
            case CG:
                timerCG.stop();
                tCG += timerCG.elapsed(TimeUnit.MILLISECONDS);
                timerCG.reset();
                break;
            case STRESS:
                timerStress.stop();
                tStress += timerStress.elapsed(TimeUnit.MILLISECONDS);
                timerStress.reset();
                break;
        }
    }

    public double getTotalTime(TimingTask task){
        switch (task){
            case BC:
                return tBC;
            case CG:
                return tCG;
            case STRESS:
                return tStress;
        }
        return  0.0;
    }

    public double getAverageTime(TimingTask task){
        switch (task){
            case BC:
                return tBC *1.0/countBC;
            case CG:
                return tCG*1.0/countCG;
            case STRESS:
                return tStress*1.0/countStress;
        }
        return  0.0;
    }

    public long[] getTotalTimeDistribution(TimingTask task)
        throws MPIException {
        LongBuffer mpiOnlyTimingBuffer =  ParallelOps.mpiOnlyBuffer;
        mpiOnlyTimingBuffer.position(0);
        long [] mpiOnlyTimingArray = new long[ParallelOps.worldProcsCount];
        switch (task){
            case BC:
                mpiOnlyTimingBuffer.put(tBC);
                break;
            case CG:
                mpiOnlyTimingBuffer.put(tCG);
                break;
            case STRESS:
                mpiOnlyTimingBuffer.put(tStress);
                break;
        }
        ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
        return mpiOnlyTimingArray;
    }

    public long[] getCountDistribution(TimingTask task) throws MPIException{
        LongBuffer mpiOnlyTimingBuffer =  ParallelOps.mpiOnlyBuffer;
        mpiOnlyTimingBuffer.position(0);
        long [] mpiOnlyTimingArray = new long[ParallelOps.worldProcsCount];
        switch (task){
            case BC:
                mpiOnlyTimingBuffer.put(countBC);
                break;
            case CG:
                mpiOnlyTimingBuffer.put(countCG);
                break;
            case STRESS:
                mpiOnlyTimingBuffer.put(countStress);
                break;
        }
        ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
        mpiOnlyTimingBuffer.position(0);
        mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
        return mpiOnlyTimingArray;
    }


}
