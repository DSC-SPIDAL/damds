package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

public class BCInternalTimings {
    public static enum TimingTask{
        BOFZ, MM
    }

    public static void init(int numThreads){
        timerBofZ = new Stopwatch[numThreads];
        IntStream.range(0, numThreads).forEach(i -> timerBofZ[i] = Stopwatch.createUnstarted());
        tBofZ = new long[numThreads];
        countBofZ = new long[numThreads];

        timerMM = new Stopwatch[numThreads];
        IntStream.range(0,numThreads).forEach(i -> timerMM[i] = Stopwatch.createUnstarted());
        tMM = new long[numThreads];
        countMM = new long[numThreads];
    }

    private static Stopwatch [] timerBofZ;
    private static Stopwatch [] timerMM;

    private static long [] tBofZ;
    private static long [] tMM;

    private static long [] countBofZ;
    private static long [] countMM;

    public static void startTiming(TimingTask task, int threadIdx){
        switch (task){
            case BOFZ:
                timerBofZ[threadIdx].start();
                ++countBofZ[threadIdx];
                break;
            case MM:
                timerMM[threadIdx].start();
                ++countMM[threadIdx];
                break;
        }
    }

    public static void endTiming(TimingTask task, int threadIdx){
        switch (task){
            case BOFZ:
                timerBofZ[threadIdx].stop();
                tBofZ[threadIdx] += timerBofZ[threadIdx].elapsed(TimeUnit.MILLISECONDS);
                timerBofZ[threadIdx].reset();
                break;
            case MM:
                timerMM[threadIdx].stop();
                tMM[threadIdx] += timerMM[threadIdx].elapsed(
                    TimeUnit.MILLISECONDS);
                timerMM[threadIdx].reset();
                break;
        }
    }

    public static double getTotalTime(TimingTask task){
        switch (task){
            case BOFZ:
                return Arrays.stream(tBofZ).reduce(0, (i,j) -> i+j);
            case MM:
                return Arrays.stream(tMM).reduce(0, (i,j) -> i+j);
        }
        return  0.0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case BOFZ:
                return Arrays.stream(tBofZ).reduce(0, (i,j) -> i+j) *1.0 / Arrays.stream(countBofZ).reduce(0, (i,j)->i+j);
            case MM:
                return Arrays.stream(tMM).reduce(0, (i,j) -> i+j) *1.0 / Arrays.stream(countMM).reduce(0, (i,j)->i+j);
        }
        return  0.0;
    }


}
