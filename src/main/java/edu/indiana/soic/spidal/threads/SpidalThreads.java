package edu.indiana.soic.spidal.threads;

import edu.indiana.soic.spidal.damds.Utils;

import java.util.*;
import java.util.concurrent.*;

public class SpidalThreads {
    private List<AbstractWorker> workers;
    private ExecutorService executor;
    private Map<Integer, Queue<TaskFuture>> queues;
    private Map<Integer, BlockingQueue<TaskFuture>> blockingQueues;
    private boolean busy;
    private int threadCount;
    private boolean bind;
    private int totalCores;
    private int startingCore;

    public SpidalThreads(int threadCount, boolean busy, boolean bind, int totalCores, int startingCore) {
        this.busy = busy;
        this.threadCount = threadCount;
        this.bind = bind;
        this.totalCores = totalCores;
        this.startingCore = startingCore;

        this.executor = Executors.newFixedThreadPool(threadCount);
        this.queues = new HashMap<>();
        this.blockingQueues = new HashMap<>();
        this.workers = new ArrayList<>();
        // start the workers
        startWorkers(threadCount, busy);
    }

    private void startWorkers(int threadCount, boolean busy) {
        for (int i = 0; i < threadCount; i++) {
            AbstractWorker worker;
            if (busy) {
                Queue<TaskFuture> queue = new LinkedList<>();
                queues.put(i, queue);
                worker = new BusyWorker(i, queue);
            } else {
                BlockingQueue<TaskFuture> queue = new ArrayBlockingQueue<>(8);
                blockingQueues.put(i, queue);
                worker = new BlockingQueueWorker(i, queue);
            }
            if (this.bind) {
                worker.setBind(bind, totalCores, startingCore);
            }
            workers.add(worker);
            executor.submit(worker);
        }
    }

    public void shutDown() {
        for (AbstractWorker w : workers) {
            w.stop();
        }
        executor.shutdown();
    }

    public void execTask(int threadIdx, Task<Integer> task) {
        if (busy) {
            submitTaskToBusyQueue(task, null, threadIdx);
        } else {
            submitTaskToBlockingQueue(task, null, threadIdx);
        }
    }

    /**
     * Submits the same code executed by many threads, changing based on task id
     */
    public void execTasks(MultiTask multiTask) {
        multiTask.run();
    }

    public void forall(Task<Integer> task) {
        CountDownLatch latch = new CountDownLatch(threadCount);

        if (this.busy) {
            for (int i = 0; i < threadCount; i++) {
                submitTaskToBusyQueue(task, latch, i);
            }
        } else {
            for (int i = 0; i < threadCount; i++) {
                submitTaskToBlockingQueue(task, latch, i);
            }
        }

        try {
            latch.await();
        } catch (InterruptedException e) {
            Utils.printAndThrowRuntimeException(new RuntimeException(e));
        }
    }

    private void submitTaskToBlockingQueue(Task<Integer> task, CountDownLatch latch, int i) {
        BlockingQueue<TaskFuture> queue = blockingQueues.get(i);
        TaskFuture future = new TaskFuture(i, task, latch);
        try {
            queue.put(future);
        } catch (InterruptedException e) {
            Utils.printAndThrowRuntimeException(new RuntimeException(e));
        }
    }

    private void submitTaskToBusyQueue(Task<Integer> task, CountDownLatch latch, int i) {
        Queue<TaskFuture> queue = queues.get(i);
        TaskFuture future = new TaskFuture(i, task, latch);
        queue.add(future);
    }
}
