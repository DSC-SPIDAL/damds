package edu.indiana.soic.spidal.damds.threads;

import java.util.Arrays;
import java.util.BitSet;

public class ThreadBitAssigner {
    /* Hard coded values for Juliet*/
    private static int spn = 2; // sockets per node
    private static int htpc = 2; // hyper threads per core

    public static void main(String[] args) {
        /* The 24 core patterns are
         * 1x24
         * 2x12
         * 3x8
         * 4x6
         * 6x4
         * 8x3
         * 12x2
         * 24x1
         */
        /*int cps = 12;*/
        /*int[] tppArray = new int[]{1,2,3,4,6,8,12,24};*/

        /* The 36 core patterns are
         * 1x36
         * 2x18
         * 3x12
         * 4x9
         * 6x6
         * 9x4
         * 12x3
         * 18x2
         * 36x1
         */
        int cps = 18;
        int[] tppArray = new int[]{1,2,3,4,6,9,12,18,36};
        int nodes = 1;

        int cpn = cps * spn;
        for (int tpp : tppArray) {
            int ppn = cpn / tpp; // process per node
            int procs = ppn * nodes;
            System.out.println("-----" + tpp + "x" + ppn + "x" + nodes + "-----");
            for (int rank = 0; rank < procs; ++rank) {
                System.out.println("  Rank: " + rank);
                for (int threadIdx = 0; threadIdx < tpp; ++threadIdx) {
                    int[] bitset = getBitMask(rank, threadIdx, tpp, cps);
                    System.out.println("    Thread: " + threadIdx + "  " +Arrays.toString(bitset));
                }
            }
        }
    }

    private static int[] getBitMask(int rank, int threadIdx, int tpp, int cps){
        int cpn = cps * spn;
        int ppn = cpn / tpp; // process per node

        // Assuming continuous ranking within a node
        int nodeLocalRank = rank % ppn;

        int[] bitset = new int[htpc];
        int idx = 0;
        for (int j = 0; j < htpc; ++j) {
            bitset[idx++] = nodeLocalRank * tpp + threadIdx + (cpn * j);
        }

        return bitset;
    }

    public static BitSet getBitSet(int rank, int threadIdx, int tpp, int cps){
        int cpn = cps * spn;
        int[] bitMask = getBitMask(rank, threadIdx, tpp, cps);
        BitSet bitSet = new BitSet(cpn);
        for(int mask: bitMask){
            bitSet.set(mask);
        }
        return bitSet;
    }
}
