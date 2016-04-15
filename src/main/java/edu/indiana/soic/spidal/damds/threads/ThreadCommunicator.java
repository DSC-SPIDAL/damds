package edu.indiana.soic.spidal.damds.threads;

import edu.indiana.soic.spidal.common.DoubleStatistics;
import edu.indiana.soic.spidal.common.RefObj;
import edu.indiana.soic.spidal.damds.ParallelOps;
import net.openhft.lang.io.Bytes;

import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;

public class ThreadCommunicator {
    private int threadCount;
    private int[] intBuffer;
    private double[] doubleBuffer;
    private double[] pointsBuffer;
    private DoubleStatistics[] doubleStatisticsBuffer;
    private CyclicBarrier barrier;

    public ThreadCommunicator(int threadCount, int numberDataPoints, int targetDimension) {
        this.threadCount = threadCount;
        intBuffer = new int[threadCount];
        doubleBuffer = new double[threadCount];
        pointsBuffer = new double[numberDataPoints*targetDimension];
        doubleStatisticsBuffer = new DoubleStatistics[threadCount];
        for (int i = 0; i < threadCount; ++i){
            doubleStatisticsBuffer[i] = new DoubleStatistics();
        }
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
    public void sumIntOverThreads(int threadIdx, RefObj<Integer> val)
        throws BrokenBarrierException, InterruptedException {
        intBuffer[threadIdx] = val.getValue();
        barrier.await();
        int sum = 0;
        for (int i = 0; i < threadCount; ++i){
            sum += intBuffer[i];
        }
        val.setValue(sum);
    }

    public void sumDoublesOverThreads(int threadIdx, RefObj<Double> val)
        throws BrokenBarrierException, InterruptedException {
        doubleBuffer[threadIdx] = val.getValue();
        barrier.await();
        double sum = 0.0;
        for (int i = 0; i < threadCount; ++i){
            sum += doubleBuffer[i];
        }
        val.setValue(sum);
    }

    public void sumDoubleStatisticsOverThreads(int threadIdx, DoubleStatistics val)
        throws BrokenBarrierException, InterruptedException {
        doubleStatisticsBuffer[threadIdx].copyFrom(val);
        barrier.await();
        DoubleStatistics sum = doubleStatisticsBuffer[0];
        if (threadIdx == 0){
            for (int i = 1; i < threadCount; ++i){
                sum.combine(doubleStatisticsBuffer[i]);
            }
        }
        barrier.await();
        val.copyFrom(sum);
    }

    public void bcastIntOverThreads(int threadIdx, RefObj<Integer> val, int root)
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

    public void bcastDoubleOverThreads(int threadIdx, RefObj<Double> val, int root)
        throws BrokenBarrierException, InterruptedException {
        if (threadIdx == root){
            final double value = val.getValue();
            for (int i = 0; i < threadCount; ++i){
                doubleBuffer[i] = value;
            }
        }
        barrier.await();
        val.setValue(doubleBuffer[threadIdx]);
    }

    public void bcastDoubleStatisticsOverThreads(int threadIdx, DoubleStatistics val, int root)
        throws BrokenBarrierException, InterruptedException {
        if (threadIdx == root){
            for (int i = 0; i < threadCount; ++i){
                doubleStatisticsBuffer[i].copyFrom(val);
            }
        }
        barrier.await();
        val.copyFrom(doubleStatisticsBuffer[threadIdx]);
    }

    public void barrier()
        throws BrokenBarrierException, InterruptedException {
        barrier.await();
    }


    public void bcastDoubleArrayOverThreads(int threadIdx, double[] preX, int root)
        throws BrokenBarrierException, InterruptedException {
        if (threadIdx == root){
            System.arraycopy(preX, 0, pointsBuffer, 0, preX.length);
        }
        barrier.await();
        System.arraycopy(pointsBuffer, 0, preX, 0, pointsBuffer.length);
    }

    public synchronized void collect(
        int startIndex, double[] val, Bytes bytes) {
        bytes.position(startIndex);
        for (double aVal : val) {
            bytes.writeDouble(aVal);
        }
    }

    public synchronized void copy(Bytes from, double[] to, int count) {
        from.position(0);
        for (int i = 0; i < count; ++i){
            to[i] = from.readDouble();
        }

    }

    public synchronized void collect(int startIndex, double[] val, Bytes
            bytes, int threadId) {
        int pos = startIndex;
        for (int i =0; i < val.length; ++i) {
            bytes.position(pos);
            bytes.writeDouble(val[i]);
            if (ParallelOps.worldProcRank == 1 && threadId == 1) {
                if (pos == ((8013 - 5000)*3+2)*Double.BYTES){
                    System.out.println("************* val=" + val[i] + " " +
                            "buffer=" + bytes.readDouble(pos));
                }
            }
            pos+=Double.BYTES;
        }

        if (ParallelOps.worldProcRank == 1 && threadId == 1) {
            System.out.println("######## startIdx=" +startIndex + " val" +
                    ".length=" + val.length + " last pos=" + (pos-Double
                    .BYTES));
        }

        /*if (ParallelOps.worldProcRank == 0 && threadId == 1) {
            System.out.println("++Rank=" + ParallelOps.worldProcRank + " " +
                    "Tid=" +
                    "" + threadId + " inCollect startIdx=" + startIndex +
                    " mmapXWriteBytes[2600][1]=" + ParallelOps
                    .mmapXWriteBytes.readDouble((2600 * 3 + 1)*Double.BYTES)
                    + " val[2600][1]=" + val[(2600-2500)*3+1]);
        }

        if (ParallelOps.worldProcRank == 1 && threadId == 0) {
            System.out.println("++Rank=" + ParallelOps.worldProcRank + " " +
                    "Tid=" +
                    "" + threadId + " inCollect startIdx=" + startIndex +
                    " mmapXWriteBytes[7200][2]=" + ParallelOps
                    .mmapXWriteBytes.readDouble(((7200 - 5000) * 3 + 2))
                    *Double.BYTES + " val[7200][2]=" + val[(7200-5000)*3+2]);
        }*/

        if (ParallelOps.worldProcRank == 1 && threadId == 1) {
            System.out.println("++Rank=" + ParallelOps.worldProcRank + " " +
                    "Tid=" +
                    "" + threadId + " inCollect startIdx=" + startIndex +
                    " mmapXWriteBytes[8013][2]=" + ParallelOps
                    .mmapXWriteBytes.readDouble(((8013 - 5000) * 3 + 2))
                    *Double.BYTES + " val[8013][2]=" + val[(8013-7500)*3+2]);
        }
    }
}
