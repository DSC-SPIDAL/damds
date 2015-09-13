package edu.indiana.soic.spidal.damds;

import edu.indiana.soic.spidal.common.DoubleStatistics;
import edu.indiana.soic.spidal.common.Range;
import edu.indiana.soic.spidal.common.RangePartitioner;
import mpi.*;
import net.openhft.lang.io.ByteBufferBytes;
import net.openhft.lang.io.Bytes;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;
import java.nio.channels.FileChannel;
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
    public static int mmapLeadCgProcRank;
    public static int mmapLeadCgProcCount;
    public static int mmapLeadWorldProcRank;
    public static int mmapLeadWorldProcRankLocalToNode;
    public static int mmapProcsRowCount;

    // mmap leaders form one communicating group and the others (followers)
    // belong to another communicating group.
    public static Intracomm cgProcComm;
    public static int[] mmapLeadsXRowCounts;
    public static int[] mmapLeadsXByteExtents;
    public static int[] mmapLeadsXDisplas;

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
    static DoubleBuffer pointBuffer;
    private static ByteBuffer statBuffer;
    private static DoubleBuffer doubleBuffer;
    private static IntBuffer intBuffer;
    private static IntBuffer twoIntBuffer;
    public static LongBuffer threadsAndMPIBuffer;
    public static LongBuffer mpiOnlyBuffer;

    public static Bytes mmapXWriteBytes;
    public static Bytes mmapXReadBytes;
    public static ByteBuffer mmapXReadByteBuffer;
    public static Bytes fullXBytes;
    public static ByteBuffer fullXByteBuffer;
    public static Bytes[] fullXBytesSlices;
    public static ByteBuffer[] fullXByteBufferSlices;
    public static Win fullXByteBufferWindow;

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
        mmapLeadWorldProcRankLocalToNode =
            isMmapLead
                ? worldProcRankLocalToNode
                : (q * mmapIdLocalToNode + (mmapIdLocalToNode < r
                                                ? mmapIdLocalToNode
                                                : r));
        mmapLeadWorldProcRank = worldProcRank - (worldProcRankLocalToNode
                                             - mmapLeadWorldProcRankLocalToNode);
        for (int i = 0; i < mmapProcsCount; ++i){
            mmapProcsWorldRanks[i] = mmapLeadWorldProcRank +i;
        }

        // Create mmap leaders' communicator
        cgProcComm = worldProcsComm.split(isMmapLead ? 0 : 1, worldProcRank);
        if (!isMmapLead){
            cgProcComm = null;
        }

        // Communicator for processes within a  memory map group
        mmapProcComm = worldProcsComm.split((nodeId*mmapsPerNode)+mmapIdLocalToNode, worldProcRank);

        /* Allocate basic buffers for communication */
        statBuffer = MPI.newByteBuffer(DoubleStatistics.extent);
        doubleBuffer = MPI.newDoubleBuffer(1);
        intBuffer = MPI.newIntBuffer(1);
        twoIntBuffer = MPI.newIntBuffer(2);

        machineName = MPI.getProcessorName();
        parallelPattern =
            "---------------------------------------------------------\n"
            + "Machine:" + machineName + ' ' + threadCount + 'x'
            + worldProcsPerNode + 'x' + nodeCount;
        Utils.printMessage(parallelPattern);
    }

    public static void tearDownParallelism() throws MPIException {
        fullXByteBufferWindow.free();
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

        // TODO - should be able to remove this using fullXByteBuffer
        // Allocate a point buffer
        pointBuffer = MPI.newDoubleBuffer(globalRowCount * targetDimension);

        // Allocate timing buffers
        mpiOnlyBuffer = MPI.newLongBuffer(worldProcsCount);
        threadsAndMPIBuffer = MPI.newLongBuffer(worldProcsCount * threadCount);

        twoIntBuffer.put(0, isMmapLead ? cgProcComm.getRank() : -1);
        twoIntBuffer.put(1, isMmapLead ? cgProcComm.getSize() : -1);
        mmapProcComm.bcast(twoIntBuffer, 2, MPI.INT, 0);
        mmapLeadCgProcRank = twoIntBuffer.get(0);
        mmapLeadCgProcCount = twoIntBuffer.get(1);


        mmapProcsRowCount = IntStream.range(mmapLeadWorldProcRank,
                                            mmapLeadWorldProcRank + mmapProcsCount)
            .map(i -> procRowRanges[i].getLength())
            .sum();
        mmapLeadsXRowCounts = new int[mmapLeadCgProcCount];
        mmapLeadsXByteExtents = new int[mmapLeadCgProcCount];
        mmapLeadsXDisplas = new int[mmapLeadCgProcCount];
        if (isMmapLead){
            mmapLeadsXRowCounts[mmapLeadCgProcRank] = mmapProcsRowCount;
            cgProcComm.allGather(mmapLeadsXRowCounts, 1, MPI.INT);
            for (int i = 0; i < mmapLeadCgProcCount; ++i){
                mmapLeadsXByteExtents[i] = mmapLeadsXRowCounts[i] * targetDimension * Double.BYTES;
            }


        }
        mmapProcComm.bcast(mmapLeadsXByteExtents, mmapLeadCgProcCount, MPI.INT, 0);
        mmapLeadsXDisplas[0] = 0;
        System.arraycopy(mmapLeadsXByteExtents, 0, mmapLeadsXDisplas, 1,
                         mmapLeadCgProcCount - 1);
        Arrays.parallelPrefix(mmapLeadsXDisplas, (m, n) -> m + n);

        final String fullXFname = machineName + ".mmapId." + mmapIdLocalToNode +".fullX.bin";
        try (FileChannel fullXFc = FileChannel.open(Paths.get(mmapScratchDir,
                                                             fullXFname),
                                                   StandardOpenOption.CREATE,StandardOpenOption.WRITE,StandardOpenOption.READ)) {

            int mmapXWriteByteExtent = procRowCount * targetDimension * Double.BYTES;
            long mmapXWriteByteOffset = (procRowStartOffset - procRowRanges[mmapLeadWorldProcRank].getStartIndex()) * targetDimension * Double.BYTES;
            int fullXByteExtent = globalRowCount * targetDimension * Double.BYTES;
            long fullXByteOffset = 0L;

            fullXBytes = ByteBufferBytes.wrap(fullXFc.map(
                FileChannel.MapMode.READ_WRITE, fullXByteOffset,
                fullXByteExtent));
            fullXByteBuffer = fullXBytes.sliceAsByteBuffer(fullXByteBuffer);
            fullXByteBufferWindow = new Win(fullXByteBuffer, fullXByteExtent, Byte.BYTES, MPI.INFO_NULL, cgProcComm);
            fullXBytesSlices = new Bytes[mmapLeadCgProcCount];
            fullXByteBufferSlices = new ByteBuffer[mmapLeadCgProcCount];
            for (int i = 0; i < mmapLeadCgProcCount; ++i){
                final int offset = mmapLeadsXDisplas[i];
                int length = mmapLeadsXByteExtents[i];
                fullXBytesSlices[i] = fullXBytes.slice(offset, length);
                fullXByteBufferSlices[i] = fullXBytesSlices[i].sliceAsByteBuffer(fullXByteBufferSlices[i]);
            }
            mmapXReadBytes = fullXBytesSlices[mmapLeadCgProcRank];
            mmapXReadByteBuffer = mmapXReadBytes.sliceAsByteBuffer(mmapXReadByteBuffer);
            mmapXWriteBytes = mmapXReadBytes.slice(mmapXWriteByteOffset,
                                                    mmapXWriteByteExtent);
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
        cgProcComm.allGatherv(mmapXReadByteBuffer,
                              mmapLeadsXByteExtents[mmapLeadCgProcRank], MPI.BYTE,
                              fullXByteBuffer, mmapLeadsXByteExtents,
                              mmapLeadsXDisplas, MPI.BYTE);
    }

    public static void iPartialXAllGather() throws MPIException{
        Request req = cgProcComm.iAllGatherv(mmapXReadByteBuffer,
                                       mmapLeadsXByteExtents[mmapLeadCgProcRank],
                                       MPI.BYTE, fullXByteBuffer,
                                       mmapLeadsXByteExtents, mmapLeadsXDisplas,
                                       MPI.BYTE);
        req.waitFor();
        req.free();
    }

    public static void partialXAllGatherLinearRing() throws MPIException {
        final int mmapLeadsSub1 = mmapLeadCgProcCount - 1;
        int recvFromRank = (mmapLeadCgProcRank + mmapLeadsSub1) % mmapLeadCgProcCount;
        int sendToRank = (mmapLeadCgProcRank+1)% mmapLeadCgProcCount;
        int recvFullXSliceIdx = recvFromRank;
        int sendFullXSliceIdx = mmapLeadCgProcRank;
        for (int i = 0; i < mmapLeadsSub1; ++i){
            cgProcComm.sendRecv(fullXByteBufferSlices[sendFullXSliceIdx],
                                mmapLeadsXByteExtents[sendFullXSliceIdx],
                                MPI.BYTE, sendToRank, sendToRank,
                                fullXByteBufferSlices[recvFullXSliceIdx],
                                mmapLeadsXByteExtents[recvFullXSliceIdx],
                                MPI.BYTE, recvFromRank, mmapLeadCgProcRank);
            sendFullXSliceIdx = recvFullXSliceIdx;
            recvFullXSliceIdx = (recvFullXSliceIdx + mmapLeadsSub1) % mmapLeadCgProcCount;
        }
    }

    public static void partialXAllGatherRemoteMemory() throws MPIException {

        fullXByteBufferWindow.fence(0);
        for (int i = 0; i < mmapLeadCgProcCount; ++i) {
            if (i == mmapLeadCgProcRank) continue;
            fullXByteBufferWindow.put(mmapXReadByteBuffer,
                                      mmapLeadsXByteExtents[mmapLeadCgProcRank],
                                      MPI.BYTE, i,
                                      mmapLeadsXDisplas[mmapLeadCgProcRank],
                                      mmapLeadsXByteExtents[mmapLeadCgProcRank],
                                      MPI.BYTE);
        }
        fullXByteBufferWindow.fence(0);
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
