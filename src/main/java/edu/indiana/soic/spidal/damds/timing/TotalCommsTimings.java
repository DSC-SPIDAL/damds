package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;

import java.util.concurrent.TimeUnit;

public class TotalCommsTimings {
    public static enum TimingTask {
        COMM, ALL
    }

    private Stopwatch timerComm = Stopwatch.createUnstarted();
    private Stopwatch timerAll = Stopwatch.createUnstarted();
    private long tComp;
    private long tAll;
    private long countComp;
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
        }
    }

    public double getTotalTime(TimingTask task) {
        switch (task) {
            case COMM:
                return tComp;
            case ALL:
                return tAll;
        }
        return 0.0;
    }

    public double getAverageTime(TimingTask task) {
        switch (task) {
            case COMM:
                return tComp * 1.0 / countComp;
            case ALL:
                return tAll * 1.0 / countAll;
        }
        return 0.0;
    }

}
