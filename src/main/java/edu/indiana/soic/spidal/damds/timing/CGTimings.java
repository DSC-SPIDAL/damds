package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;

import java.util.concurrent.TimeUnit;

public class CGTimings {
    public static enum TimingTask{
        MM, INNER_PROD, CG_LOOP

    }

    private Stopwatch timerMM = Stopwatch.createUnstarted();
    private Stopwatch timerInnerProd = Stopwatch.createUnstarted();
    private Stopwatch timerCGLoop = Stopwatch.createUnstarted();

    private long tMM;
    private long tInnerProd;
    private long tCGLoop;

    private long countMM;
    private long countInnerProd;
    private long countCGLoop;

    public void startTiming(TimingTask task){
        switch (task){
            case MM:
                timerMM.start();
                ++countMM;
                break;
            case INNER_PROD:
                timerInnerProd.start();
                ++countInnerProd;
                break;
            case CG_LOOP:
                timerCGLoop.start();
                ++countCGLoop;
                break;
        }
    }

    public void endTiming(TimingTask task){
        switch (task){
            case MM:
                timerMM.stop();
                tMM += timerMM.elapsed(TimeUnit.MILLISECONDS);
                timerMM.reset();
                break;
            case INNER_PROD:
                timerInnerProd.stop();
                tInnerProd += timerInnerProd.elapsed(TimeUnit.MILLISECONDS);
                timerInnerProd.reset();
                break;
            case CG_LOOP:
                timerCGLoop.stop();
                tCGLoop += timerCGLoop.elapsed(TimeUnit.MILLISECONDS);
                timerCGLoop.reset();
                break;
        }
    }

    public double getTotalTime(TimingTask task){
        switch (task){
            case MM:
                return tMM;
            case INNER_PROD:
                return tInnerProd;
            case CG_LOOP:
                return tCGLoop;
        }
        return  0.0;
    }

    public double getAverageTime(TimingTask task){
        switch (task){
            case MM:
                return tMM *1.0/ countMM;
            case INNER_PROD:
                return tInnerProd *1.0/ countInnerProd;
            case CG_LOOP:
                return tCGLoop *1.0/ countCGLoop;
        }
        return  0.0;
    }



}
