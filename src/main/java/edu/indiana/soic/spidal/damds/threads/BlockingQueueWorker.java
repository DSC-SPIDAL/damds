package edu.indiana.soic.spidal.damds.threads;

import edu.indiana.soic.spidal.damds.Utils;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.CountDownLatch;

public class BlockingQueueWorker extends AbstractWorker {
    private BlockingQueue<TaskFuture> blockingQueue;

    public BlockingQueueWorker(int threadIdx, BlockingQueue<TaskFuture> queue) {
        super(threadIdx);
        this.blockingQueue = queue;
    }

    @Override
    public void run() {
        try {
            while (this.run) {
                TaskFuture t = blockingQueue.take();
                if (t != null) {
                    Task<Integer> task = t.getTask();
                    task.run(t.getThreadId());

                    CountDownLatch latch = t.getLatch();
                    if (latch != null) {
                        latch.countDown();
                    }
                }
            }
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }
}
