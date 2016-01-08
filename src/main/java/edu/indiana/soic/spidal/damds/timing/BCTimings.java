package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.damds.ParallelOps;
import mpi.MPIException;

import java.nio.LongBuffer;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

public class BCTimings {
    public static enum TimingTask{
        BC_INTERNAL,COMM, BC_MERGE, BC_EXTRACT
    }

    private static int numThreads;
    public static void init(int numThreads){
        timerBCInternal = new Stopwatch[numThreads];
        IntStream.range(0,numThreads).forEach(i -> timerBCInternal[i] = Stopwatch.createUnstarted());
        tBCInternal = new long[numThreads];
        countBCInternal = new long[numThreads];
        BCTimings.numThreads = numThreads;
    }

    private static Stopwatch [] timerBCInternal;
    private static Stopwatch timerComm = Stopwatch.createUnstarted();
    private static Stopwatch timerBCMerge = Stopwatch.createUnstarted();
    private static Stopwatch timerBCExtract = Stopwatch.createUnstarted();

    private static long [] tBCInternal;
    private static long tComm;
    private static long tBCMerge;
    private static long tBCExtract;

    private static long [] countBCInternal;
    private static long countComm;
    private static long countBCMerge;
    private static long countBCExtract;

    public static void startTiming(TimingTask task, int threadIdx){
        switch (task){
            case BC_INTERNAL:
                timerBCInternal[threadIdx].start();
                ++countBCInternal[threadIdx];
                break;
            case COMM:
                timerComm.start();
                ++countComm;
                break;
            case BC_MERGE:
                timerBCMerge.start();
                ++countBCMerge;
                break;
            case BC_EXTRACT:
                timerBCExtract.start();
                ++countBCExtract;
                break;
        }
    }

    public static void endTiming(TimingTask task, int threadIdx){
        switch (task){
            case BC_INTERNAL:
                timerBCInternal[threadIdx].stop();
                tBCInternal[threadIdx] += timerBCInternal[threadIdx].elapsed(TimeUnit.MILLISECONDS);
                timerBCInternal[threadIdx].reset();
                break;
            case COMM:
                timerComm.stop();
                tComm += timerComm.elapsed(TimeUnit.MILLISECONDS);
                timerComm.reset();
                break;
            case BC_MERGE:
                timerBCMerge.stop();
                tBCMerge += timerBCMerge.elapsed(TimeUnit.MILLISECONDS);
                timerBCMerge.reset();
                break;
            case BC_EXTRACT:
                timerBCExtract.stop();
                tBCExtract += timerBCExtract.elapsed(TimeUnit.MILLISECONDS);
                timerBCExtract.reset();
                break;
        }
    }

    public static double getTotalTime(TimingTask task){
        switch (task){
            case BC_INTERNAL:
                return Arrays.stream(tBCInternal).reduce(0, (i,j) -> i+j);
            case COMM:
                return tComm;
            case BC_MERGE:
                return tBCMerge;
            case BC_EXTRACT:
                return tBCExtract;
        }
        return  0.0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case BC_INTERNAL:
                return Arrays.stream(tBCInternal).reduce(0, (i,j) -> i+j) *1.0 / Arrays.stream(countBCInternal).reduce(0, (i,j)->i+j);
            case COMM:
                return tComm *1.0/ countComm;
            case BC_MERGE:
                return tBCMerge * 1.0 / countBCMerge;
            case BC_EXTRACT:
                return tBCExtract * 1.0 / countBCExtract;
        }
        return  0.0;
    }

    public static long[] getTotalTimeDistribution(TimingTask task)
        throws MPIException {
        LongBuffer threadsAndMPITimingBuffer =
            ParallelOps.threadsAndMPIBuffer;
        LongBuffer mpiOnlyTimingBuffer =  ParallelOps.mpiOnlyBuffer;
        threadsAndMPITimingBuffer.position(0);
        mpiOnlyTimingBuffer.position(0);
        long [] threadsAndMPITimingArray = new long[numThreads * ParallelOps.worldProcsCount];
        long [] mpiOnlyTimingArray = new long[ParallelOps.worldProcsCount];
        switch (task){
            case BC_INTERNAL:
                threadsAndMPITimingBuffer.put(tBCInternal);
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
            case BC_MERGE:
                mpiOnlyTimingBuffer.put(tBCMerge);
                ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
                mpiOnlyTimingBuffer.position(0);
                mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
                return mpiOnlyTimingArray;
            case BC_EXTRACT:
                mpiOnlyTimingBuffer.put(tBCExtract);
                ParallelOps.gather(mpiOnlyTimingBuffer, 1, 0);
                mpiOnlyTimingBuffer.position(0);
                mpiOnlyTimingBuffer.get(mpiOnlyTimingArray);
                return mpiOnlyTimingArray;
        }
        return null;
    }


}
