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

    private Stopwatch timerComp = Stopwatch.createUnstarted();
    private long tComp;
    private long countComp;

    public void startTiming(TimingTask task, int threadIdx){
        switch (task){
            case COMP:
                timerComp.start();
                ++countComp;
                break;
        }
    }

    public void endTiming(TimingTask task, int threadIdx){
        switch (task){
            case COMP:
                timerComp.stop();
                tComp += timerComp.elapsed(TimeUnit.MILLISECONDS);
                timerComp.reset();
                break;
        }
    }

    public double getTotalTime(TimingTask task){
        switch (task){
            case COMP:
                return tComp;
        }
        return  0.0;
    }

    public double getAverageTime(TimingTask task){
        switch (task){
            case COMP:
                return tComp * 1.0/ countComp;
        }
        return  0.0;
    }

    /*public long[] getTotalTimeDistribution(TimingTask task)
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
    }*/


}
