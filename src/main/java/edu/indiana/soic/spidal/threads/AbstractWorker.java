package edu.indiana.soic.spidal.threads;

import net.openhft.affinity.Affinity;

import java.util.BitSet;

public abstract class AbstractWorker implements Runnable{
    protected int threadIdx;

    protected boolean bind;

    protected int totalCores;

    protected int core;

    protected boolean run = true;

    public AbstractWorker(int threadIdx) {
        this.threadIdx = threadIdx;
    }

    public void setBind(boolean bind, int totalCores, int core) {
        this.bind = bind;
        this.totalCores = totalCores;
        this.core = core;
    }

    public static void bind(int totalCores, int core) {
        BitSet bitSet = new BitSet(totalCores);
        bitSet.set(core);
        Affinity.setAffinity(bitSet);
    }

    public void stop() {
        this.run = false;
    }
}
