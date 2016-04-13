package edu.indiana.soic.spidal.damds.threads;

import edu.indiana.soic.spidal.common.DoubleStatistics;
import edu.indiana.soic.spidal.common.RefObj;

import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.CyclicBarrier;

public class ThreadCommunicator {
    private int threadCount;
    private int[] intBuffer;
    private DoubleStatistics[] doubleStatisticsBuffer;
    private CyclicBarrier barrier;

    public ThreadCommunicator(int threadCount) {
        this.threadCount = threadCount;
        intBuffer = new int[threadCount];
        doubleStatisticsBuffer = new DoubleStatistics[threadCount];
        barrier = new CyclicBarrier(threadCount);
    }

    /**
     * Sum value over threads and collects at root
     * @param threadIdx the thread index
     * @param val the value
     * @return the summation
     * @throws BrokenBarrierException
     * @throws InterruptedException
     */
    public void sumOverThreads(int threadIdx, RefObj<Integer> val)
        throws BrokenBarrierException, InterruptedException {
        intBuffer[threadIdx] = val.getValue();
        barrier.await();
        int sum = 0;
        for (int i = 0; i < threadCount; ++i){
            sum += intBuffer[i];
        }
        val.setValue(sum);
    }

    public void sumOverThreads(int threadIdx, DoubleStatistics val)
        throws BrokenBarrierException, InterruptedException {
        doubleStatisticsBuffer[threadIdx].copyFrom(val);
        barrier.await();
        DoubleStatistics sum = doubleStatisticsBuffer[0];
        if (threadIdx == 0){
            for (int i = 0; i < threadCount; ++i){
                sum.combine(doubleStatisticsBuffer[i]);
            }
        }
        barrier.await();
        val.copyFrom(sum);
    }

    public void bcastOverThreads(int threadIdx, RefObj<Integer> val, int root)
        throws BrokenBarrierException, InterruptedException {
        if (threadIdx == root){
            final Integer value = val.getValue();
            for (int i = 0; i < threadCount; ++i){
                intBuffer[i] = value;
            }
        }
        barrier.await();
        val.setValue(intBuffer[threadIdx]);
    }

    public void bcastOverThreads(int threadIdx, DoubleStatistics val, int root)
        throws BrokenBarrierException, InterruptedException {
        if (threadIdx == root){
            for (int i = 0; i < threadCount; ++i){
                doubleStatisticsBuffer[i].copyFrom(val);
            }
        }
        barrier.await();
        val.copyFrom(doubleStatisticsBuffer[threadIdx]);
    }

    public void threadBarrier()
        throws BrokenBarrierException, InterruptedException {
        barrier.await();
    }


}
