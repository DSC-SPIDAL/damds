package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.damds.ParallelOps;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.Arrays;
import java.util.OptionalLong;
import java.util.concurrent.TimeUnit;

public class BCInternalTimings {
    public static enum TimingTask {
        BOFZ, MM
    }

    private Stopwatch timerBofZ;
    private Stopwatch timerMM;

    private long tBofZ;
    private long tMM;

    private long countBofZ;
    private long countMM;

    public void startTiming(TimingTask task) {
        switch (task) {
            case BOFZ:
                timerBofZ.start();
                ++countBofZ;
                break;
            case MM:
                timerMM.start();
                ++countMM;
                break;
        }
    }

    public void endTiming(TimingTask task) {
        switch (task) {
            case BOFZ:
                timerBofZ.stop();
                tBofZ +=
                    timerBofZ.elapsed(TimeUnit.MILLISECONDS);
                timerBofZ.reset();
                break;
            case MM:
                timerMM.stop();
                tMM +=
                    timerMM.elapsed(TimeUnit.MILLISECONDS);
                timerMM.reset();
                break;
        }
    }

    public double getTotalTime(TimingTask task) {
        switch (task) {
            case BOFZ:
                return  tBofZ;
            case MM:
                return tMM;
        }
        return 0.0;
    }

    public double getAverageTime(TimingTask task) {
        switch (task) {
            case BOFZ:
                return tBofZ * 1.0 / countBofZ;
            case MM:
                return tMM * 1.0 / countMM;
        }
        return 0.0;
    }

    /*public long[] getTotalTimeDistribution(TimingTask task)
        throws MPIException {
        LongBuffer threadsAndMPITimingBuffer = ParallelOps.threadsAndMPIBuffer;
        threadsAndMPITimingBuffer.position(0);
        long[] array = new long[numThreads * ParallelOps.worldProcsCount];
        switch (task) {
            case BOFZ:
                threadsAndMPITimingBuffer.put(tBofZ);
                break;
            case MM:
                threadsAndMPITimingBuffer.put(tMM);
                break;
        }
        ParallelOps.gather(threadsAndMPITimingBuffer, numThreads, 0);
        threadsAndMPITimingBuffer.position(0);
        threadsAndMPITimingBuffer.get(array);
        return array;
    }*/


}
