package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;

import java.util.concurrent.TimeUnit;

public class TotalCommsTimings {
    public static enum TimingTask {
        COMM, ALL, STATS, STRESS, BARRIER
    }

    private Stopwatch timerComm = Stopwatch.createUnstarted();
    private Stopwatch timerStress = Stopwatch.createUnstarted();
    private Stopwatch timerStats = Stopwatch.createUnstarted();
    private Stopwatch timerBar = Stopwatch.createUnstarted();
    private Stopwatch timerAll = Stopwatch.createUnstarted();
    private long tComp;
    private long tStress;
    private long tStats;
    private long tBar;
    private long tAll;
    private long countComp;
    private long countStress;
    private long countStats;
    private long countBar;
    private long countAll;

    public void startTiming(TimingTask task) {
        switch (task) {
            case COMM:
                timerComm.start();
                ++countComp;
                break;
            case ALL:
                timerAll.start();
                ++countAll;
                break;
            case STATS:
                timerStats.start();
                ++countStats;
                break;
            case STRESS:
                timerStress.start();
                ++countStress;
                break;
            case BARRIER:
                timerBar.start();
                ++countBar;
                break;
        }
    }

    public void endTiming(TimingTask task) {
        switch (task) {
            case COMM:
                timerComm.stop();
                tComp += timerComm.elapsed(TimeUnit.MILLISECONDS);
                timerComm.reset();
                break;
            case ALL:
                timerAll.stop();
                tAll += timerAll.elapsed(TimeUnit.MILLISECONDS);
                timerAll.reset();
                break;
            case STRESS:
                timerStress.stop();
                tStress += timerStress.elapsed(TimeUnit.MILLISECONDS);
                timerStress.reset();
                break;
            case STATS:
                timerStats.stop();
                tStats += timerStats.elapsed(TimeUnit.MILLISECONDS);
                timerStats.reset();
                break;
            case BARRIER:
                timerBar.stop();
                tBar += timerBar.elapsed(TimeUnit.MILLISECONDS);
                timerBar.reset();
                break;
        }
    }

    public double getTotalTime(TimingTask task) {
        switch (task) {
            case COMM:
                return tComp;
            case ALL:
                return tAll;
            case STATS:
                return tStats;
            case STRESS:
                return tStress;
            case BARRIER:
                return tBar;
        }
        return 0.0;
    }

    public double getAverageTime(TimingTask task) {
        switch (task) {
            case COMM:
                return tComp * 1.0 / countComp;
            case ALL:
                return tAll * 1.0 / countAll;
            case STRESS:
                return tStress * 1.0 / countStress;
            case STATS:
                return tStats * 1.0 / countStats;
            case BARRIER:
                return tBar * 1.0 / countBar;
        }
        return 0.0;
    }

}
