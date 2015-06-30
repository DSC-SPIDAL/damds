package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

public class MMTimings {
    public static enum TimingTask{
        MM_INTERNAL,COMM
    }

    public static void init(int numThreads){
        timerMMInternal = new Stopwatch[numThreads];
        IntStream.range(0, numThreads).forEach(i -> timerMMInternal[i] = Stopwatch.createUnstarted());
        tMMInternal = new long[numThreads];
        countMMInternal = new long[numThreads];
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
                return Arrays.stream(tMMInternal).reduce(0, (i,j) -> i+j);
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

}
