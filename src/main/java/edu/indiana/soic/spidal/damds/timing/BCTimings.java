package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.damds.ParallelOps;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

public class BCTimings {
    public static enum TimingTask{
        BC_INTERNAL,COMM, BC_MERGE, BC_EXTRACT
    }

    private Stopwatch timerBCInternal = Stopwatch.createUnstarted();
    private Stopwatch timerComm = Stopwatch.createUnstarted();
    private Stopwatch timerBCMerge = Stopwatch.createUnstarted();
    private Stopwatch timerBCExtract = Stopwatch.createUnstarted();

    private long tBCInternal;
    private long tComm;
    private long tBCMerge;
    private long tBCExtract;

    private long countBCInternal;
    private long countComm;
    private long countBCMerge;
    private long countBCExtract;

    public void startTiming(TimingTask task){
        switch (task){
            case BC_INTERNAL:
                timerBCInternal.start();
                ++countBCInternal;
                break;
            case COMM:
                timerComm.start();
                ++countComm;
                break;
            case BC_MERGE:
                timerBCMerge.start();
                ++countBCMerge;
                break;
            case BC_EXTRACT:
                timerBCExtract.start();
                ++countBCExtract;
                break;
        }
    }

    public void endTiming(TimingTask task){
        switch (task){
            case BC_INTERNAL:
                timerBCInternal.stop();
                tBCInternal += timerBCInternal.elapsed(TimeUnit.MILLISECONDS);
                timerBCInternal.reset();
                break;
            case COMM:
                timerComm.stop();
                tComm += timerComm.elapsed(TimeUnit.MILLISECONDS);
                timerComm.reset();
                break;
            case BC_MERGE:
                timerBCMerge.stop();
                tBCMerge += timerBCMerge.elapsed(TimeUnit.MILLISECONDS);
                timerBCMerge.reset();
                break;
            case BC_EXTRACT:
                timerBCExtract.stop();
                tBCExtract += timerBCExtract.elapsed(TimeUnit.MILLISECONDS);
                timerBCExtract.reset();
                break;
        }
    }

    public double getTotalTime(TimingTask task){
        switch (task){
            case BC_INTERNAL:
                return tBCInternal;
            case COMM:
                return tComm;
            case BC_MERGE:
                return tBCMerge;
            case BC_EXTRACT:
                return tBCExtract;
        }
        return  0.0;
    }

    public double getAverageTime(TimingTask task){
        switch (task){
            case BC_INTERNAL:
                return tBCInternal * 1.0 / countBCInternal;
            case COMM:
                return tComm *1.0/ countComm;
            case BC_MERGE:
                return tBCMerge * 1.0 / countBCMerge;
            case BC_EXTRACT:
                return tBCExtract * 1.0 / countBCExtract;
        }
        return  0.0;
    }

    /*public long[] getTotalTimeDistribution(TimingTask task)
        throws MPIException {
        LongBuffer threadsAndMPITimingBuffer =
            ParallelOps.threadsAndMPIBuffer;
        LongBuffer mpiOnlyTimingBuffer =  ParallelOps.mpiOnlyBuffer;
        threadsAndMPITimingBuffer.position(0);
        mpiOnlyTimingBuffer.position(0);
        long [] threadsAndMPITimingArray = new long[numThreads * ParallelOps.worldProcsCount];
        long [] mpiOnlyTimingArray = new long[ParallelOps.worldProcsCount];
        switch (task){
            case BC_INTERNAL:
                threadsAndMPITimingBuffer.put(tBCInternal);
                ParallelOps.gather(threadsAndMPITimingBuffer, numThreads, 0);
                threadsAndMPITimingBuffer.position(0);
                threadsAndMPITimingBuffer.get(threadsAndMPITimingArray);
                return threadsAndMPITimingArray;
            case COMM:
                mpiOnlyTimingBuffer.put(tComm);
                ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
                mpiOnlyTimingBuffer.position(0);
                mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
                return mpiOnlyTimingArray;
            case BC_MERGE:
                mpiOnlyTimingBuffer.put(tBCMerge);
                ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
                mpiOnlyTimingBuffer.position(0);
                mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
                return mpiOnlyTimingArray;
            case BC_EXTRACT:
                mpiOnlyTimingBuffer.put(tBCExtract);
                ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
                mpiOnlyTimingBuffer.position(0);
                mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
                return mpiOnlyTimingArray;
        }
        return null;
    }*/
}
