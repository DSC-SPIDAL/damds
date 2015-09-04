package edu.indiana.soic.spidal.damds;

import edu.indiana.soic.spidal.common.DoubleStatistics;
import edu.indiana.soic.spidal.common.Range;
import edu.indiana.soic.spidal.common.RangePartitioner;
import mpi.Intracomm;
import mpi.MPI;
import mpi.MPIException;
import net.openhft.lang.io.ByteBufferBytes;
import net.openhft.lang.io.Bytes;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.stream.IntStream;

public class ParallelOps {
    public static String machineName;
    public static int nodeCount=1;
    public static int threadCount=1;

    public static int nodeId;

    public static Intracomm worldProcsComm;
    public static int worldProcRank;
    public static int worldProcsCount;
    public static int worldProcsPerNode;

    public static Intracomm mmapProcComm;
    // Number of memory mapped groups per process
    public static int mmapsPerNode;
    public static String mmapScratchDir;
    public static int worldProcRankLocalToNode;
    public static int mmapIdLocalToNode;
    public static int mmapProcsCount;
    public static boolean isMmapLead;
    public static int[] mmapProcsWorldRanks;
    public static int mmapLeadWorldRank;
    public static int mmapLeadWorldRankLocalToNode;

    // mmap leaders form one communicating group and the others (followers)
    // belong to another communicating group.
    public static Intracomm cgProcComm;
    public static int cgProcRank;
    public static int cgProcsCount;
    public static int[] cgProcsRowCounts;
    public static int[] cgProcsPartialXByteExtents;
    public static int[] cgProcsPartialXDisplas;

    public static String parallelPattern;
    public static Range[] procRowRanges;
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

    public static Bytes partialXLeaderReadBytes;
    public static ByteBuffer partialXLeaderReadByteBuffer;
    public static Bytes partialXWriteBytes;
    public static Bytes fullXBytes;
    public static ByteBuffer fullXByteBuffer;
 /*   public static Bytes lockAndCountBytes;
    public static final int LOCK_OFFSET = 0;
    public static final int COUNT_OFFSET = 4;
    private static final int LOCK_AND_COUNT_EXTENT = 12; // 12 bytes = 4 byte lock + 8 byte (long) count*/


    public static void setupParallelism(String[] args) throws MPIException {
        MPI.Init(args);
        worldProcsComm = MPI.COMM_WORLD; //initializing MPI world communicator
        worldProcRank = worldProcsComm.getRank();
        worldProcsCount = worldProcsComm.getSize();



        /* Create communicating groups */
        worldProcsPerNode = worldProcsCount / nodeCount;
        if ((worldProcsPerNode * nodeCount) != worldProcsCount) {
            Utils.printAndThrowRuntimeException(
                "Inconsistent MPI counts Nodes " + nodeCount + " Size "
                + worldProcsCount);
        }

        worldProcRankLocalToNode = worldProcRank % worldProcsPerNode;
        nodeId = worldProcRank / worldProcsPerNode;
        int q = worldProcsPerNode / mmapsPerNode;
        int r = worldProcsPerNode % mmapsPerNode;

        // Memory mapped groups and communicating groups
        mmapIdLocalToNode =
            worldProcRankLocalToNode < r * (q + 1)
                ? worldProcRankLocalToNode / (q + 1)
                : (worldProcRankLocalToNode - r) / q;
        mmapProcsCount = worldProcRankLocalToNode < r*(q+1) ? q+1 : q;
        isMmapLead = worldProcRankLocalToNode % mmapProcsCount == 0;
        mmapProcsWorldRanks = new int[mmapProcsCount];
        mmapLeadWorldRankLocalToNode =
            isMmapLead
                ? worldProcRankLocalToNode
                : (q * mmapIdLocalToNode + (mmapIdLocalToNode < r
                                                ? mmapIdLocalToNode
                                                : r));
        mmapLeadWorldRank = worldProcRank - (worldProcRankLocalToNode
                                             - mmapLeadWorldRankLocalToNode);
        for (int i = 0; i < mmapProcsCount; ++i){
            mmapProcsWorldRanks[i] = mmapLeadWorldRank +i;
        }

        // Leaders talk, their color is 0
        // Followers will get a communicator of color 1,
        // but will make sure they don't talk ha ha :)
        cgProcComm = worldProcsComm.split(isMmapLead ? 0 : 1, worldProcRank);
        cgProcRank = cgProcComm.getRank();
        cgProcsCount = cgProcComm.getSize();

        // Communicator for processes within a  memory map group
        mmapProcComm = worldProcsComm.split((nodeId*mmapsPerNode)+mmapIdLocalToNode, worldProcRank);

        /* Allocate basic buffers for communication */
        statBuffer = MPI.newByteBuffer(DoubleStatistics.extent);
        doubleBuffer = MPI.newDoubleBuffer(1);
        intBuffer = MPI.newIntBuffer(1);

        machineName = MPI.getProcessorName();
        parallelPattern =
            "---------------------------------------------------------\n"
            + "Machine:" + machineName + ' ' + threadCount + 'x'
            + worldProcsPerNode + 'x' + nodeCount;
        Utils.printMessage(parallelPattern);
    }

    public static void tearDownParallelism() throws MPIException {
        // End MPI
        MPI.Finalize();
    }

    public static void setParallelDecomposition(int globalRowCount, int targetDimension)
        throws IOException, MPIException {
        //	First divide points among processes
        procRowRanges = RangePartitioner.partition(globalRowCount,
                                                       worldProcsCount);
        Range rowRange = procRowRanges[worldProcRank]; // The range of points for this process

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
        IntStream.range(0, threadCount)
            .parallel()
            .forEach(threadIdx -> {
                         Range threadRowRange = threadRowRanges[threadIdx];
                         threadRowCounts[threadIdx] =
                             threadRowRange.getLength();
                         threadRowStartOffsets[threadIdx] =
                             threadRowRange.getStartIndex();
                         threadPointStartOffsets[threadIdx] =
                             threadRowStartOffsets[threadIdx] * globalColCount;
                     });

        // Allocate vector buffers
        partialPointBuffer = MPI.newDoubleBuffer(procRowCount * targetDimension);
        pointBuffer = MPI.newDoubleBuffer(globalRowCount * targetDimension);
        mpiOnlyBuffer = MPI.newLongBuffer(worldProcsCount);
        threadsAndMPIBuffer = MPI.newLongBuffer(worldProcsCount * threadCount);

        cgProcsRowCounts = new int[cgProcsCount];
        cgProcsPartialXByteExtents = new int[cgProcsCount];
        cgProcsPartialXDisplas = new int[cgProcsCount];
        if (isMmapLead){
            int rowCount = IntStream.range(mmapLeadWorldRank,
                                           mmapLeadWorldRank + mmapProcsCount)
                .map(i -> procRowRanges[i].getLength())
                .sum();
            cgProcsRowCounts[cgProcRank] = rowCount;
            cgProcComm.allGather(cgProcsRowCounts, 1, MPI.INT);
            for (int i = 0; i < cgProcsCount; ++i){
                cgProcsPartialXByteExtents[i] = cgProcsRowCounts[i] * targetDimension * Double.BYTES;
            }

            cgProcsPartialXDisplas[0] = 0;
            System.arraycopy(cgProcsPartialXByteExtents, 0, cgProcsPartialXDisplas, 1, cgProcsCount - 1);
            Arrays.parallelPrefix(cgProcsPartialXDisplas, (m, n) -> m + n);
        }


        final String partialXFname = machineName + ".mmapId." + mmapIdLocalToNode + ".partialX.bin";
        final String fullXFname = machineName + ".mmapId." + mmapIdLocalToNode +".fullX.bin";
        final String lockAndCountFname = machineName + ".lockAndCount.bin";
        try (FileChannel partialXFc = FileChannel.open(Paths.get(mmapScratchDir,
                                                                 partialXFname),
                                                       StandardOpenOption
                                                           .CREATE,
                                                       StandardOpenOption.READ,
                                                       StandardOpenOption
                                                           .WRITE);
            FileChannel fullXFc = FileChannel.open(Paths.get(mmapScratchDir,
                                                             fullXFname),
                                                   StandardOpenOption.CREATE,StandardOpenOption.WRITE,StandardOpenOption.READ);
            /*FileChannel lockAndCountFc = FileChannel.open(Paths.get(
                                                              mmapScratchDir,
                                                              lockAndCountFname),
                                                          StandardOpenOption
                                                              .CREATE,
                                                          StandardOpenOption.READ,
                                                          StandardOpenOption.WRITE)*/){


            long partialXLeaderReadExtent = cgProcsRowCounts[cgProcRank] * targetDimension * Double.BYTES;
            long partialXWriteExtent = procRowCount * targetDimension * Double.BYTES;
            long partialXOffset = (procRowStartOffset - procRowRanges[mmapLeadWorldRank].getStartIndex()) * targetDimension * Double.BYTES;
            long fullXExtent = globalRowCount * targetDimension * Double.BYTES;
            long fullXOffset = 0L;

            partialXWriteBytes = ByteBufferBytes.wrap(
                partialXFc.map(FileChannel.MapMode.READ_WRITE, partialXOffset,
                               partialXWriteExtent));

            partialXLeaderReadBytes = isMmapLead ? ByteBufferBytes.wrap(partialXFc.map(
                FileChannel.MapMode.READ_ONLY, partialXOffset, partialXLeaderReadExtent)) : null;
            partialXLeaderReadByteBuffer = isMmapLead ? partialXLeaderReadBytes.sliceAsByteBuffer(partialXLeaderReadByteBuffer) : null;


            fullXBytes = ByteBufferBytes.wrap(fullXFc.map(FileChannel.MapMode.READ_WRITE,
                                            fullXOffset, fullXExtent));
            fullXByteBuffer = fullXBytes.sliceAsByteBuffer(fullXByteBuffer);

            /*lockAndCountBytes = ByteBufferBytes.wrap(lockAndCountFc.map(
                FileChannel.MapMode.READ_WRITE, 0, LOCK_AND_COUNT_EXTENT));*/

            // Print debug info in order of world ranks
            for (int i = 0; i < worldProcsCount; ++i){
                intBuffer.put(0, i);
                worldProcsComm.bcast(intBuffer, 1, MPI.INT, 0);
                int next = intBuffer.get(0);
                if (next == worldProcRank){
                    try (BufferedWriter bw = Files.newBufferedWriter(Paths.get("mmap.debug" + worldProcRank + ".out.txt"), StandardOpenOption.CREATE, StandardOpenOption.WRITE)) {

                        PrintWriter writer = new PrintWriter(bw, true);
                        // Good it's my turn to print
                        writer.println(
                            "World rank: " + worldProcRank + " on " + machineName);
                        writer.println("  mmapIdLocalToNode:             " + mmapIdLocalToNode);
                        writer.println("  mmapProcsCount:                " + mmapProcsCount);
                        writer.println("  isMmapLead:                    " + isMmapLead);
                        writer.println("  mmapProcsWorldRanks:           " + Arrays.toString(
                            mmapProcsWorldRanks));
                        writer.println("  mmapLeadWorldRankLocalToNode:  "
                                           + "" + mmapLeadWorldRankLocalToNode);
                        writer.println("  mmapLeadWorldRank:             " + mmapLeadWorldRank);
                        writer.println("  cgProcRank:                    " + cgProcRank);
                        writer.println("  cgProcsCount:                  "
                                           + "" + cgProcsCount);
                        writer.println("  cgProcsRowCounts:              "
                                       + Arrays.toString(cgProcsRowCounts));
                        writer.println("  partialXLeaderReadExtent:      "
                                       + partialXLeaderReadExtent);
                        writer.println("  partialXWriteExtent:           "
                                       + partialXWriteExtent);
                        writer.println("  partialXOffset:                "
                                       + partialXOffset);
                        writer.println("  fullXExtent:                   "
                                       + fullXExtent);
                        writer.println("  fullXOffset:                   "
                                       + fullXOffset);
                    }
                }
            }
        }
    }

    public static DoubleStatistics allReduce(DoubleStatistics stat)
        throws MPIException {
        stat.addToBuffer(statBuffer, 0);
        worldProcsComm.allReduce(statBuffer, DoubleStatistics.extent, MPI.BYTE,
                                 DoubleStatistics.reduceSummaries());
        return DoubleStatistics.getFromBuffer(statBuffer, 0);
    }

    public static double allReduce(double value) throws MPIException{
        doubleBuffer.put(0, value);
        worldProcsComm.allReduce(doubleBuffer, 1, MPI.DOUBLE, MPI.SUM);
        return doubleBuffer.get(0);
    }

    public static int allReduce(int value) throws MPIException{
        intBuffer.put(0, value);
        worldProcsComm.allReduce(intBuffer, 1, MPI.INT, MPI.SUM);
        return intBuffer.get(0);
    }

    public static void partialXAllGather() throws MPIException {
        partialXLeaderReadByteBuffer.position(0);
        fullXByteBuffer.position(0);
        cgProcComm.allGatherv(partialXLeaderReadByteBuffer,
                              cgProcsPartialXByteExtents[cgProcRank],
                              MPI.BYTE, fullXByteBuffer,
                              cgProcsPartialXByteExtents,
                              cgProcsPartialXDisplas, MPI.BYTE);
//        fullXBytes.force();
    }

    public static void broadcast(DoubleBuffer buffer, int extent, int root)
        throws MPIException {
        worldProcsComm.bcast(buffer, extent, MPI.DOUBLE, root);
    }

    public static void gather(LongBuffer buffer, int count, int root)
        throws MPIException {
        worldProcsComm.gather(buffer, count, MPI.LONG, root);
    }
}
