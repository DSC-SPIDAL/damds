package edu.indiana.soic.spidal.threads;

import java.util.Queue;
import java.util.concurrent.CountDownLatch;

public class BusyWorker extends AbstractWorker {
    private Queue<TaskFuture> taskQueue;

    public BusyWorker(int threadIdx, Queue<TaskFuture> queue) {
        super(threadIdx);
        this.taskQueue = queue;
    }

    @Override
    public void run() {
        if (bind) {
            bind(totalCores, core);
        }

        while (true) {
            TaskFuture t = taskQueue.poll();
            if (t != null) {
                System.out.println("Exec");
                Task<Integer> task = t.getTask();
                task.run(t.getThreadId());

                CountDownLatch latch = t.getLatch();
                if (latch != null) {
                    latch.countDown();
                }
            }
        }
    }
}
