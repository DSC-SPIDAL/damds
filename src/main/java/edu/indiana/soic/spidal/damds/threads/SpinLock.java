package edu.indiana.soic.spidal.damds.threads;

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
}