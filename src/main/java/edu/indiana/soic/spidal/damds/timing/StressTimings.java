package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.damds.ParallelOps;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

public class StressTimings {
    public enum TimingTask{
        STRESS_INTERNAL,COMM, STRESS_MERGE, STRESS_EXTRACT
    }

    private Stopwatch timerStressInternal = Stopwatch.createUnstarted();
    private Stopwatch timerComm = Stopwatch.createUnstarted();
    private Stopwatch timerBCMerge = Stopwatch.createUnstarted();
    private Stopwatch timerBCExtract = Stopwatch.createUnstarted();

    private long tStressInternal;
    private long tComm;
    private long tBCMerge;
    private long tBCExtract;

    private long countStressInternal;
    private long countComm;
    private long countBCMerge;
    private long countBCExtract;

    public void startTiming(TimingTask task, int threadIdx){
        switch (task){
            case STRESS_INTERNAL:
                timerStressInternal.start();
                ++countStressInternal;
                break;
            case COMM:
                timerComm.start();
                ++countComm;
                break;
            case STRESS_MERGE:
                timerBCMerge.start();
                ++countBCMerge;
                break;
            case STRESS_EXTRACT:
                timerBCExtract.start();
                ++countBCExtract;
                break;
        }
    }

    public void endTiming(TimingTask task, int threadIdx){
        switch (task){
            case STRESS_INTERNAL:
                timerStressInternal.stop();
                tStressInternal += timerStressInternal.elapsed(TimeUnit.MILLISECONDS);
                timerStressInternal.reset();
                break;
            case COMM:
                timerComm.stop();
                tComm += timerComm.elapsed(TimeUnit.MILLISECONDS);
                timerComm.reset();
                break;
            case STRESS_MERGE:
                timerBCMerge.stop();
                tBCMerge += timerBCMerge.elapsed(TimeUnit.MILLISECONDS);
                timerBCMerge.reset();
                break;
            case STRESS_EXTRACT:
                timerBCExtract.stop();
                tBCExtract += timerBCExtract.elapsed(TimeUnit.MILLISECONDS);
                timerBCExtract.reset();
                break;
        }
    }

    public double getTotalTime(TimingTask task){
        switch (task){
            case STRESS_INTERNAL:
                return tStressInternal;
            case COMM:
                return tComm;
            case STRESS_MERGE:
                return tBCMerge;
            case STRESS_EXTRACT:
                return tBCExtract;
        }
        return  0.0;
    }

    public double getAverageTime(TimingTask task){
        switch (task){
            case STRESS_INTERNAL:
                return tStressInternal * 1.0 / countStressInternal;
            case COMM:
                return tComm *1.0/ countComm;
            case STRESS_MERGE:
                return tBCMerge * 1.0 / countBCMerge;
            case STRESS_EXTRACT:
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
            case STRESS_INTERNAL:
                threadsAndMPITimingBuffer.put(tStressInternal);
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
            case STRESS_MERGE:
                mpiOnlyTimingBuffer.put(tBCMerge);
                ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
                mpiOnlyTimingBuffer.position(0);
                mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
                return mpiOnlyTimingArray;
            case STRESS_EXTRACT:
                mpiOnlyTimingBuffer.put(tBCExtract);
                ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
                mpiOnlyTimingBuffer.position(0);
                mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
                return mpiOnlyTimingArray;
        }
        return null;
    }*/


}
