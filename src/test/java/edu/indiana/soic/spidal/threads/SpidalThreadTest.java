package edu.indiana.soic.spidal.threads;

import edu.indiana.soic.spidal.damds.ParallelOps;
import junit.framework.TestCase;

public class SpidalThreadTest extends TestCase {

    public void testA() {
        SpidalThreads threads = new SpidalThreads(10, false, false, 48, 1);
        threads.forall((threadId) -> print(threadId));

        threads.forall((threadId) -> print(threadId));
        System.out.println("Finished");
    }

    public void print(int thread) {
        for (int i = 0; i < 10000; i++) {
            for (int j = 0; j < 1000; j++) {
                int k = i / (j + 1);
            }
        }
        System.out.println(thread);
    }
}
