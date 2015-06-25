package edu.indiana.soic.spidal.damds;

public class Utils {
    public static void printAndThrowRuntimeException(RuntimeException e) {
        e.printStackTrace(System.out);
        throw e;
    }

    public static void printAndThrowRuntimeException(String message) {
        System.out.println(message);
        throw new RuntimeException(message);
    }

    public static void printMessage(String msg) {
        if (ParallelOps.procRank != 0) {
            return;
        }
        System.out.println(msg);
    }
}
