package edu.indiana.soic.spidal.damds;

public class Utils {

    int threadId;

    public Utils(int threadId) {
        this.threadId = threadId;
    }

    public void printAndThrowRuntimeException(RuntimeException e) {
        e.printStackTrace(System.out);
        throw e;
    }

    public void printAndThrowRuntimeException(String message) {
        System.out.println(message);
        throw new RuntimeException(message);
    }

    public void printMessage(String msg) {
        if (ParallelOps.worldProcRank != 0 || threadId != 0) {
            return;
        }
        System.out.println(msg);
    }
}
