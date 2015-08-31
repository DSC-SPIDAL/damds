package edu.indiana.soic.spidal.damds;

import edu.indiana.soic.spidal.common.DoubleStatistics;
import edu.indiana.soic.spidal.common.Range;
import edu.indiana.soic.spidal.common.RangePartitioner;
import mpi.Intracomm;
import mpi.MPI;
import mpi.MPIException;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.stream.IntStream;

public class ParallelOps {
    public static String machineName;
    public static int nodeCount=1;
    public static int threadCount=1;

    public static Intracomm worldProcComm;
    public static int worldProcRank;
    public static int worldProcCount;

    public static int worldProcsPerNode;
    // Number of communicating groups per process
    public static int cgPerNode;
    public static String cgScratchDir;
    public static int worldNodeLocalRank;
    public static int cgIdLocalToNode;
    public static int cgProcCountPerNode;
    public static boolean isCgLead;
    public static boolean isCgIOLead;
    public static Intracomm cgComm;
    public static int cgProcRank;
    public static int cgProcCount;
    public static int[] cgPeersWorldRanks;
    public static int cgLeadWorldRank;
    public static int cgLeadNodeLocalRank;


    public static String parallelPattern;
    public static Range procRowRange;
    public static int procRowStartOffset;
    public static int procRowCount;

    public static long procPointStartOffset;
    public static Range[] threadRowRanges;
    public static int[] threadRowStartOffsets;
    public static int[] threadRowCounts;


    public static int[] threadPointStartOffsets;

    public static int globalColCount;

    // Buffers for MPI operations
    private static ByteBuffer statBuffer;
    private static DoubleBuffer doubleBuffer;
    private static IntBuffer intBuffer;
    static DoubleBuffer partialPointBuffer;
    static DoubleBuffer pointBuffer;
    public static LongBuffer threadsAndMPIBuffer;
    public static LongBuffer mpiOnlyBuffer;

    public static FileChannel cgPartialPoints;


    public static void setupParallelism(String[] args) throws MPIException {
        MPI.Init(args);
        worldProcComm = MPI.COMM_WORLD; //initializing MPI world communicator
        worldProcRank = worldProcComm.getRank();
        worldProcCount = worldProcComm.getSize();

        worldProcsPerNode = worldProcCount / nodeCount;

        if ((worldProcsPerNode * nodeCount) != worldProcCount) {
            Utils.printAndThrowRuntimeException(
                "Inconsistent MPI counts Nodes " + nodeCount + " Size " +
                worldProcCount);
        }

        /* Create communicating groups */
        worldProcsPerNode = worldProcCount / nodeCount;

        worldNodeLocalRank = worldProcRank % worldProcsPerNode;
        int q = worldProcsPerNode / cgPerNode;
        int r = worldProcsPerNode % cgPerNode;

        // Communicating group
        cgIdLocalToNode = worldNodeLocalRank < r*(q+1) ? worldNodeLocalRank/(q+1) : (worldNodeLocalRank-r)/q;
        cgProcCountPerNode = worldNodeLocalRank < r*(q+1) ? q+1 : q;
        isCgLead = worldNodeLocalRank % cgProcCountPerNode == 0;
        isCgIOLead = worldNodeLocalRank == 0;
        cgPeersWorldRanks = new int[cgProcCountPerNode];
        cgLeadNodeLocalRank = isCgLead ? worldNodeLocalRank : (q*cgIdLocalToNode + (cgIdLocalToNode < r ? cgIdLocalToNode : r));
        cgLeadWorldRank = worldProcRank - (worldNodeLocalRank - cgLeadNodeLocalRank);
        for (int i = 0; i < cgProcCountPerNode; ++i){
            cgPeersWorldRanks[i] = cgLeadWorldRank+i;
        }


        // Leaders talk, their color is 0
        // Followers will get a communicator of color 1, but will make sure they don't talk ha ha :)
        cgComm = worldProcComm.split(isCgLead ? 0 : 1, worldProcRank);
        cgProcRank = cgComm.getRank();
        cgProcCount = cgComm.getSize();

        /* Allocate basic buffers for communication */
        statBuffer = MPI.newByteBuffer(DoubleStatistics.extent);
        doubleBuffer = MPI.newDoubleBuffer(1);
        intBuffer = MPI.newIntBuffer(1);

        machineName = MPI.getProcessorName();
        parallelPattern =
            "---------------------------------------------------------\n" +
            "Machine:" + machineName + ' ' +
            threadCount + 'x' + worldProcsPerNode + 'x' + nodeCount;
        Utils.printMessage(parallelPattern);
    }

    public static void tearDownParallelism() throws MPIException {
        // End MPI
        MPI.Finalize();
    }

    public static void setParallelDecomposition(int globalRowCount, int targetDimension) {
        //	First divide points among processes
        Range[] rowRanges = RangePartitioner.partition(globalRowCount,
                                                       worldProcCount);
        Range rowRange = rowRanges[worldProcRank]; // The range of points for this process

        procRowRange = rowRange;
        procRowStartOffset = rowRange.getStartIndex();
        procRowCount = rowRange.getLength();
        globalColCount = globalRowCount;
        procPointStartOffset = ((long)procRowStartOffset) * globalColCount;

        // Next partition points per process among threads
        threadRowRanges = RangePartitioner.partition(procRowCount, threadCount);
        threadRowCounts = new int[threadCount];
        threadRowStartOffsets = new int[threadCount];
        threadPointStartOffsets = new int[threadCount];
        IntStream.range(0, threadCount).parallel().forEach(
            threadIdx -> {
                Range threadRowRange = threadRowRanges[threadIdx];
                threadRowCounts[threadIdx] = threadRowRange.getLength();
                threadRowStartOffsets[threadIdx] =
                    threadRowRange.getStartIndex();
                threadPointStartOffsets[threadIdx] =
                    threadRowStartOffsets[threadIdx] * globalColCount;
            });

        // Allocate vector buffers
        partialPointBuffer = MPI.newDoubleBuffer(procRowCount * targetDimension);
        pointBuffer = MPI.newDoubleBuffer(globalRowCount * targetDimension);
        mpiOnlyBuffer = MPI.newLongBuffer(worldProcCount);
        threadsAndMPIBuffer = MPI.newLongBuffer(worldProcCount * threadCount);

        // TODO - continue from here
//        cgPartialPoints = FileChannel.open(Paths.get(cgScratchDir, ))
    }

    public static DoubleStatistics allReduce(DoubleStatistics stat) throws
        MPIException {
        stat.addToBuffer(statBuffer,0);
        worldProcComm.allReduce(
            statBuffer, DoubleStatistics.extent, MPI.BYTE,
            DoubleStatistics.reduceSummaries());
        return DoubleStatistics.getFromBuffer(statBuffer, 0);
    }

    public static double allReduce(double value) throws MPIException{
        doubleBuffer.put(0, value);
        worldProcComm.allReduce(doubleBuffer, 1, MPI.DOUBLE, MPI.SUM);
        return doubleBuffer.get(0);
    }

    public static int allReduce(int value) throws MPIException{
        intBuffer.put(0, value);
        worldProcComm.allReduce(intBuffer, 1, MPI.INT, MPI.SUM);
        return intBuffer.get(0);
    }

    public static DoubleBuffer allGather(
        DoubleBuffer partialPointBuffer, int dimension) throws MPIException {

        int [] lengths = new int[worldProcCount];
        int length = procRowCount * dimension;
        lengths[worldProcRank] = length;
        worldProcComm.allGather(lengths, 1, MPI.INT);
        int [] displas = new int[worldProcCount];
        displas[0] = 0;
        System.arraycopy(lengths, 0, displas, 1, worldProcCount - 1);
        Arrays.parallelPrefix(displas, (m, n) -> m + n);
        worldProcComm.allGatherv(
            partialPointBuffer, length, MPI.DOUBLE, pointBuffer, lengths,
            displas, MPI.DOUBLE);
        return  pointBuffer;
    }

    public static void broadcast(DoubleBuffer buffer, int extent, int root)
        throws MPIException {
        worldProcComm.bcast(buffer, extent, MPI.DOUBLE, root);
    }

    public static void gather(LongBuffer buffer, int count, int root)
        throws MPIException {
        worldProcComm.gather(buffer, count, MPI.LONG, root);
    }
}
