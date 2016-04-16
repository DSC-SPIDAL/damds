package edu.indiana.soic.spidal.damds.threads;

import java.util.concurrent.CountDownLatch;

public class TaskFuture {
    private Task<Integer> task;

    private CountDownLatch latch;

    private int threadId;

    public TaskFuture(int threadId, Task<Integer> task, CountDownLatch latch) {
        this.task = task;
        this.latch = latch;
        this.threadId = threadId;
    }

    public Task<Integer> getTask() {
        return task;
    }

    public CountDownLatch getLatch() {
        return latch;
    }

    public int getThreadId() {
        return threadId;
    }
}
