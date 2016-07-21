package edu.indiana.soic.spidal.damds.threads;

import java.util.Arrays;
import java.util.BitSet;

public class ThreadBitAssigner {
    /* Hard coded values for Juliet*/
    private static int cps = 12; // cores per socket
    private static int spn = 2; // sockets per node
    private static int htpc = 2; // hyper threads per core
    private static int cpn = cps * spn; // cores per node

    public static void main(String[] args) {
        /* The patterns are
         * 1x24
         * 2x12
         * 3x8
         * 4x6
         * 6x4
         * 8x3
         * 12x2
         * 24x1
         */
        int[] tppArray = new int[]{1,2,3,4,6,8,12,24};
        int nodes = 3;

        for (int tpp : tppArray) {
            int ppn = cpn / tpp; // process per node
            int cpp = cpn / ppn; // cores per process
            int procs = ppn * nodes;
            System.out.println("-----" + tpp + "x" + ppn + "x" + nodes + "-----");
            for (int rank = 0; rank < procs; ++rank) {
                System.out.println("  Rank: " + rank);
                for (int threadIdx = 0; threadIdx < tpp; ++threadIdx) {
                    int[] bitset = getBitMask(rank, threadIdx, tpp, nodes);
                    System.out.println("    Thread: " + threadIdx + "  " +Arrays.toString(bitset));
                }
            }
        }
    }

    private static int[] getBitMask(int rank, int threadIdx, int tpp, int nodes){
        int ppn = cpn / tpp; // process per node
        int cpp = cpn / ppn; // cores per process
        int procs = ppn * nodes;

        // Assuming continuous ranking within a node
        int nodeLocalRank = rank % ppn;

        int[] bitset = new int[(cpp/tpp)*htpc];
        int idx = 0;
        for (int j = 0; j < htpc; ++j) {
            bitset[idx++] = nodeLocalRank * cpp + threadIdx + (cpn * j);
        }

        return bitset;
    }

    public static BitSet getBitSet(int rank, int threadIdx, int tpp, int nodes){
        int[] bitMask = getBitMask(rank, threadIdx, tpp, nodes);
        BitSet bitSet = new BitSet(cpn);
        for(int mask: bitMask){
            bitSet.set(mask);
        }
        return bitSet;
    }
}
