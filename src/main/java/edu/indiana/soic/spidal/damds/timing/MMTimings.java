package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.damds.ParallelOps;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.Arrays;
import java.util.OptionalLong;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

public class MMTimings {
    public enum TimingTask{
        MM_INTERNAL,COMM, MM_MERGE, MM_EXTRACT
    }

    private Stopwatch timerMMInternal = Stopwatch.createUnstarted();
    private Stopwatch timerComm = Stopwatch.createUnstarted();
    private Stopwatch timerMMMerge = Stopwatch.createUnstarted();
    private Stopwatch timerMMExtract = Stopwatch.createUnstarted();

    private long tMMInternal;
    private long tComm;
    private long tMMMerge;
    private long tMMExtract;

    private long countMMInternal;
    private long countComm;
    private long countMMMerge;
    private long countMMExtract;



    public void startTiming(TimingTask task){
        switch (task){
            case MM_INTERNAL:
                timerMMInternal.start();
                ++countMMInternal;
                break;
            case COMM:
                timerComm.start();
                ++countComm;
                break;
            case MM_MERGE:
                timerMMMerge.start();
                ++countMMMerge;
                break;
            case MM_EXTRACT:
                timerMMExtract.start();
                ++countMMExtract;
                break;
        }
    }

    public void endTiming(TimingTask task){
        switch (task){
            case MM_INTERNAL:
                timerMMInternal.stop();
                tMMInternal += timerMMInternal.elapsed(
                    TimeUnit.MILLISECONDS);
                timerMMInternal.reset();
                break;
            case COMM:
                timerComm.stop();
                tComm += timerComm.elapsed(TimeUnit.MILLISECONDS);
                timerComm.reset();
                break;
            case MM_MERGE:
                timerMMMerge.stop();
                tMMMerge += timerMMMerge.elapsed(TimeUnit.MILLISECONDS);
                timerMMMerge.reset();
                break;
            case MM_EXTRACT:
                timerMMExtract.stop();
                tMMExtract += timerMMExtract.elapsed(TimeUnit.MILLISECONDS);
                timerMMExtract.reset();
                break;
        }
    }

    public double getTotalTime(TimingTask task){
        switch (task){
            case MM_INTERNAL:
                return tMMInternal;
            case COMM:
                return tComm;
            case MM_MERGE:
                return tMMMerge;
            case MM_EXTRACT:
                return tMMExtract;
        }
        return  0.0;
    }

    public double getAverageTime(TimingTask task){
        switch (task){
            case MM_INTERNAL:
                return tMMInternal * 1.0/ countMMInternal;
            case COMM:
                return tComm *1.0/ countComm;
            case MM_MERGE:
                return tMMMerge * 1.0 / countMMMerge;
            case MM_EXTRACT:
                return tMMExtract * 1.0 / countMMExtract;
        }
        return  0.0;
    }

    /*public long[] getTotalTimeDistribution(TimingTask task)
        throws MPIException {
        LongBuffer threadsAndMPITimingBuffer =
            ParallelOps.threadsAndMPIBuffer;
        LongBuffer mpiOnlyTimingBuffer = ParallelOps.mpiOnlyBuffer;
        threadsAndMPITimingBuffer.position(0);
        mpiOnlyTimingBuffer.position(0);
        long [] threadsAndMPITimingArray = new long[numThreads * ParallelOps.worldProcsCount];
        long [] mpiOnlyTimingArray = new long[ParallelOps.worldProcsCount];
        switch (task){
            case MM_INTERNAL:
                threadsAndMPITimingBuffer.put(tMMInternal);
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
            case MM_MERGE:
                mpiOnlyTimingBuffer.put(tMMMerge);
                ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
                mpiOnlyTimingBuffer.position(0);
                mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
                return mpiOnlyTimingArray;
            case MM_EXTRACT:
                mpiOnlyTimingBuffer.put(tMMExtract);
                ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
                mpiOnlyTimingBuffer.position(0);
                mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
                return mpiOnlyTimingArray;
        }

        return null;
    }*/

    /*public long[] getCountDistribution(TimingTask task)
        throws MPIException {
        LongBuffer threadsAndMPIBuffer =
            ParallelOps.threadsAndMPIBuffer;
        LongBuffer mpiOnlyBuffer = ParallelOps.mpiOnlyBuffer;
        threadsAndMPIBuffer.position(0);
        mpiOnlyBuffer.position(0);
        long [] threadsAndMPIArray = new long[numThreads * ParallelOps.worldProcsCount];
        long [] mpiOnlyArray = new long[ParallelOps.worldProcsCount];
        switch (task){
            case MM_INTERNAL:
                threadsAndMPIBuffer.put(countMMInternal);
                ParallelOps.gather(threadsAndMPIBuffer, numThreads, 0);
                threadsAndMPIBuffer.position(0);
                threadsAndMPIBuffer.get(threadsAndMPIArray);
                return threadsAndMPIArray;
            case COMM:
                mpiOnlyBuffer.put(countComm);
                ParallelOps.gather(mpiOnlyBuffer, 1, 0);
                mpiOnlyBuffer.position(0);
                mpiOnlyBuffer.get(mpiOnlyArray);
                return mpiOnlyArray;
            case MM_MERGE:
                mpiOnlyBuffer.put(countMMMerge);
                ParallelOps.gather(mpiOnlyBuffer, 1, 0);
                mpiOnlyBuffer.position(0);
                mpiOnlyBuffer.get(mpiOnlyArray);
                return mpiOnlyArray;
            case MM_EXTRACT:
                mpiOnlyBuffer.put(countMMExtract);
                ParallelOps.gather(mpiOnlyBuffer, 1, 0);
                mpiOnlyBuffer.position(0);
                mpiOnlyBuffer.get(mpiOnlyArray);
                return mpiOnlyArray;
        }
        return null;
    }*/


}
