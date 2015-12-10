package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.damds.ParallelOps;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.Arrays;
import java.util.OptionalLong;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

public class StressInternalTimings {
    public static enum TimingTask{
        COMP
    }

    private static int numThreads;
    public static void init(int numThreads){
        timerComp = new Stopwatch[numThreads];
        IntStream.range(0, numThreads).forEach(i -> timerComp[i] = Stopwatch.createUnstarted());
        tComp = new long[numThreads];
        countComp = new long[numThreads];

        StressInternalTimings.numThreads = numThreads;
    }

    private static Stopwatch [] timerComp;

    private static long [] tComp;

    private static long [] countComp;

    public static void startTiming(TimingTask task, int threadIdx){
        switch (task){
            case COMP:
                timerComp[threadIdx].start();
                ++countComp[threadIdx];
                break;
        }
    }

    public static void endTiming(TimingTask task, int threadIdx){
        switch (task){
            case COMP:
                timerComp[threadIdx].stop();
                tComp[threadIdx] += timerComp[threadIdx].elapsed(TimeUnit.MILLISECONDS);
                timerComp[threadIdx].reset();
                break;
        }
    }

    public static double getTotalTime(TimingTask task){
        switch (task){
            case COMP:
                OptionalLong maxBofZ = Arrays.stream(tComp).max();
                return maxBofZ.isPresent() ? maxBofZ.getAsLong()*1.0 : 0.0;
        }
        return  0.0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case COMP:
                return Arrays.stream(tComp).reduce(0, (i, j) -> i + j) * 1.0 / Arrays.stream(
                    countComp).reduce(0, (i, j)-> i + j);
        }
        return  0.0;
    }

    public static long[] getTotalTimeDistribution(TimingTask task)
        throws MPIException {
        LongBuffer threadsAndMPITimingBuffer =
            ParallelOps.threadsAndMPIBuffer;
        threadsAndMPITimingBuffer.position(0);
        long [] array = new long[numThreads * ParallelOps.worldProcsCount];
        switch (task){
            case COMP:
                threadsAndMPITimingBuffer.put(tComp);
                break;
        }
        ParallelOps.gather(threadsAndMPITimingBuffer, numThreads, 0);
        threadsAndMPITimingBuffer.position(0);
        threadsAndMPITimingBuffer.get(array);
        return array;
    }


}
