package edu.indiana.soic.spidal.damds;

import java.io.*;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

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

    public static void writeOutput(double[] x, int vecLen, String outputFile)
            throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
        int N = x.length / vecLen;

        DecimalFormat format = new DecimalFormat("#.##########");
        for (int i = 0; i < N; i++) {
            int index = i * vecLen;
            writer.print(String.valueOf(i) + '\t'); // print ID.
            for (int j = 0; j < vecLen; j++) {
                writer.print(format.format(x[index + j]) + '\t'); // print
                // configuration
                // of each axis.
            }
            writer.println("1"); // print label value, which is ONE for all
            // data.
        }
        writer.flush();
        writer.close();

    }

    public static String getProcAffinityMask(int pid) throws IOException {
        byte[] bo = new byte[100];
        String pidString = String.valueOf(pid);
        String[] cmd = {"bash", "-c", "taskset -pc " + pidString};
        Process p = Runtime.getRuntime().exec(cmd);
        p.getInputStream().read(bo);
        return new String(bo);
    }
}
