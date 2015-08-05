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
        MM_INTERNAL,COMM
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

    private static long [] tMMInternal;
    private static long tComm;

    private static long [] countMMInternal;
    private static long countComm;

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
        }
    }

    public static double getTotalTime(TimingTask task){
        switch (task){
            case MM_INTERNAL:
                OptionalLong max = Arrays.stream(tMMInternal).max();
                return max.isPresent() ? max.getAsLong() * 1.0 : 0.0;
            case COMM:
                return tComm;
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
        }
        return  0.0;
    }

    public static long[] getTotalTimeDistribution(TimingTask task)
        throws MPIException {
        LongBuffer threadsAndMPITimingBuffer =
            ParallelOps.threadsAndMPITimingBuffer;
        LongBuffer mpiOnlyTimingBuffer = ParallelOps.mpiOnlyTimingBuffer;
        threadsAndMPITimingBuffer.position(0);
        mpiOnlyTimingBuffer.position(0);
        long [] threadsAndMPITimingArray = new long[numThreads * ParallelOps.procCount];
        long [] mpiOnlyTimingArray = new long[ParallelOps.procCount];
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
        }
        return null;
    }

}
