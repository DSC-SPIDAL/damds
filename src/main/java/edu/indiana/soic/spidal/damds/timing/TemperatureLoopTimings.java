package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;

import java.util.concurrent.TimeUnit;

public class TemperatureLoopTimings {
    public static enum TimingTask{
        PRE_STRESS, STRESS_LOOP
    }

    private static Stopwatch timerPreStress = Stopwatch.createUnstarted();
    private static Stopwatch timerStressLoop = Stopwatch.createUnstarted();

    private static long tPreStress;
    private static long tStressLoop;

    private static long countPreStress;
    private static long countStressLoop;

    public static void startTiming(TimingTask task){
        switch (task){
            case PRE_STRESS:
                timerPreStress.start();
                ++countPreStress;
                break;
            case STRESS_LOOP:
                timerStressLoop.start();
                ++countStressLoop;
                break;
        }
    }

    public static void endTiming(TimingTask task){
        switch (task){
            case PRE_STRESS:
                timerPreStress.stop();
                tPreStress += timerPreStress.elapsed(TimeUnit.MILLISECONDS);
                timerPreStress.reset();
                break;
            case STRESS_LOOP:
                timerStressLoop.stop();
                tStressLoop += timerStressLoop.elapsed(TimeUnit.MILLISECONDS);
                timerStressLoop.reset();
                break;
        }
    }

    public static double getTotalTime(TimingTask task){
        switch (task){
            case PRE_STRESS:
                return tPreStress;
            case STRESS_LOOP:
                return tStressLoop;
        }
        return  0.0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case PRE_STRESS:
                return tPreStress *1.0/ countPreStress;
            case STRESS_LOOP:
                return tStressLoop *1.0/ countStressLoop;
        }
        return  0.0;
    }

}
