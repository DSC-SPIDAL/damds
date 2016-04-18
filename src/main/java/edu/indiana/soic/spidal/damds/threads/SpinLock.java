package edu.indiana.soic.spidal.damds.threads;

import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.atomic.AtomicReference;

public class SpinLock
{
    private final AtomicReference<Thread> _lock = new AtomicReference<>(null);
    private final Lock _unlock = new Lock();

    public Lock lock()
    {
        Thread thread = Thread.currentThread();
        while(true)
        {
            if (!_lock.compareAndSet(null,thread))
            {
                if (_lock.get()==thread)
                    throw new IllegalStateException("SpinLock is not reentrant");
                continue;
            }
            return _unlock;
        }
    }

    public boolean isLocked()
    {
        return _lock.get()!=null;
    }

    public boolean isLockedThread()
    {
        return _lock.get()==Thread.currentThread();
    }

    public class Lock implements AutoCloseable
    {
        @Override
        public void close()
        {
            _lock.set(null);
        }
    }

    public static class Worker implements Runnable {
        private SpinLock lock;
        private CyclicBarrier cyclicBarrier;

        public Worker(SpinLock lock, CyclicBarrier cyclicBarrier) {
            this.lock = lock;
            this.cyclicBarrier = cyclicBarrier;
        }

        public Worker(SpinLock lock) {
            this.lock = lock;
        }

        @Override
        public void run() {
            Lock l = lock.lock();
            double k = 0;
                for (int i = 0; i < 1000; i++) {
                    for (int j = 0; j < 1000; j++) {
                        k = j * (k + j);
                    }
                }
            l.close();
            try {
                cyclicBarrier.await();
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (BrokenBarrierException e) {
                e.printStackTrace();
            }
            System.out.println("Finish: " + k);
        }
    }

    public static void main(String[] args) throws BrokenBarrierException, InterruptedException {
        CyclicBarrier cyclicBarrier = new CyclicBarrier(10);
        System.out.println("No waiting: " + cyclicBarrier.getNumberWaiting());
        SpinLock lock = new SpinLock();
        for (int i = 0; i < 10; i++) {
            Thread t = new Thread(new Worker(lock, cyclicBarrier));
            t.start();
        }
        Thread.sleep(10000);
        System.out.println("No waiting: " + cyclicBarrier.getNumberWaiting());
        for (int i = 0; i < 10; i++) {
            Thread t = new Thread(new Worker(lock, cyclicBarrier));
            t.start();
        }
        System.out.println("Cyclic barrier");
    }
}