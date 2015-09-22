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
    public static enum TimingTask{
        MM_INTERNAL,COMM, MM_MERGE, MM_EXTRACT
    }

    private static int numThreads;
    public static void init(int numThreads){
        timerMMInternal = new Stopwatch[numThreads];
        IntStream.range(0, numThreads).forEach(i -> timerMMInternal[i] = Stopwatch.createUnstarted());
        tMMInternal = new long[numThreads];
        countMMInternal = new long[numThreads];
        MMTimings.numThreads = numThreads;
    }

    private static Stopwatch [] timerMMInternal;
    private static Stopwatch timerComm = Stopwatch.createUnstarted();
    private static Stopwatch timerMMMerge = Stopwatch.createUnstarted();
    private static Stopwatch timerMMExtract = Stopwatch.createUnstarted();

    private static long [] tMMInternal;
    private static long tComm;
    private static long tMMMerge;
    private static long tMMExtract;

    private static long [] countMMInternal;
    private static long countComm;
    private static long countMMMerge;
    private static long countMMExtract;



    public static void startTiming(TimingTask task, int threadIdx){
        switch (task){
            case MM_INTERNAL:
                timerMMInternal[threadIdx].start();
                ++countMMInternal[threadIdx];
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

    public static void endTiming(TimingTask task, int threadIdx){
        switch (task){
            case MM_INTERNAL:
                timerMMInternal[threadIdx].stop();
                tMMInternal[threadIdx] += timerMMInternal[threadIdx].elapsed(
                    TimeUnit.MILLISECONDS);
                timerMMInternal[threadIdx].reset();
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

    public static double getTotalTime(TimingTask task){
        switch (task){
            case MM_INTERNAL:
                OptionalLong max = Arrays.stream(tMMInternal).max();
                return max.isPresent() ? max.getAsLong() * 1.0 : 0.0;
            case COMM:
                return tComm;
            case MM_MERGE:
                return tMMMerge;
            case MM_EXTRACT:
                return tMMExtract;
        }
        return  0.0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case MM_INTERNAL:
                return Arrays.stream(tMMInternal).reduce(0, (i,j) -> i+j) *1.0 / Arrays.stream(
                    countMMInternal).reduce(0, (i,j)->i+j);
            case COMM:
                return tComm *1.0/ countComm;
            case MM_MERGE:
                return tMMMerge * 1.0 / countMMMerge;
            case MM_EXTRACT:
                return tMMExtract * 1.0 / countMMExtract;
        }
        return  0.0;
    }

    public static long[] getTotalTimeDistribution(TimingTask task)
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
    }

    public static long[] getCountDistribution(TimingTask task)
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
    }


}
