package edu.indiana.soic.spidal.damds;

import edu.indiana.soic.spidal.common.DoubleStatistics;
import edu.indiana.soic.spidal.common.Range;
import edu.indiana.soic.spidal.common.RangePartitioner;
import edu.indiana.soic.spidal.damds.threads.ThreadCommunicator;
import edu.indiana.soic.spidal.damds.timing.BCTimings;
import mpi.Intracomm;
import mpi.MPI;
import mpi.MPIException;
import mpi.Op;
import net.openhft.lang.io.ByteBufferBytes;
import net.openhft.lang.io.Bytes;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.UUID;
import java.util.regex.Pattern;
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
    public static int mmapProcRank;
    public static int mmapProcsCount;
    public static boolean isMmapLead;
    public static int[] mmapProcsWorldRanks;
    public static int mmapLeadWorldRank;
    public static int mmapLeadWorldRankLocalToNode;
    public static int mmapProcsRowCount;

    // mmap leaders form one communicating group and the others (followers)
    // belong to another communicating group.
    public static Intracomm cgProcComm;
    public static int cgProcRank;
    public static int cgProcsCount;
    public static int[] cgProcsMmapRowCounts;
    public static int[] cgProcsMmapXByteExtents;
    public static int[] cgProcsMmapXDisplas;

    public static String parallelPattern;
    public static Range[] procRowRanges;
    public static Range procRowRange;
    public static int procRowStartOffset;
    public static int procRowCount;

    public static Range[] threadRowRanges;
    public static int[] threadRowStartOffsets;
    public static int[] threadRowCounts;

    public static int globalColCount;

    // Buffers for MPI operations
    private static ByteBuffer statBuffer;
    private static DoubleBuffer doubleBuffer;
    private static IntBuffer intBuffer;
    public static LongBuffer threadsAndMPIBuffer;
    public static LongBuffer mpiOnlyBuffer;

    public static String mmapLockFileNameOne;
    public static String mmapEntryLockFileName;
    public static Bytes mmapLockOne;
    public static Bytes mmapEntryLock;
    private static int FLAG = 0;
    private static int COUNT = Long.BYTES;

    public static Bytes mmapXReadBytes;
    public static ByteBuffer mmapXReadByteBuffer;
    public static Bytes mmapXWriteBytes;
    public static Bytes fullXBytes;
    public static ByteBuffer fullXByteBuffer;

    public static Bytes mmapSReadBytes;
    public static ByteBuffer mmapSReadByteBuffer;
    public static Bytes mmapSWriteBytes;

    public static ThreadCommunicator threadComm;
    
    private static Utils utils = new Utils(0);

    private static boolean isHeterogeneous;

    public static void setupParallelism(String[] args) throws MPIException {
        MPI.Init(args);
        machineName = MPI.getProcessorName();

        /* Allocate basic buffers for communication */
        statBuffer = MPI.newByteBuffer(DoubleStatistics.extent);
        doubleBuffer = MPI.newDoubleBuffer(1);
        intBuffer = MPI.newIntBuffer(1);

        worldProcsComm = MPI.COMM_WORLD; //initializing MPI world communicator
        worldProcRank = worldProcsComm.getRank();
        worldProcsCount = worldProcsComm.getSize();

        /* Create communicating groups */
        worldProcsPerNode = worldProcsCount / nodeCount;
        isHeterogeneous = (worldProcsPerNode * nodeCount) != worldProcsCount;
        if (isHeterogeneous) {
            utils.printMessage("Running in heterogeneous mode");
        }

        /* Logic to identify how many processes are within a node and
        *  the q and r values. These are used to processes to mmap groups
        *  within a node.
        *
        *  Important: the code assumes continues rank distribution
        *  within a node. */
        int[] qr = findQandR();
        int q = qr[0];
        int r = qr[1];

        // Memory mapped groups and communicating groups
        mmapIdLocalToNode =
            worldProcRankLocalToNode < r * (q + 1)
                ? worldProcRankLocalToNode / (q + 1)
                : (worldProcRankLocalToNode - r) / q;
        mmapProcsCount = worldProcRankLocalToNode < r*(q+1) ? q+1 : q;


        // Communicator for processes within a  memory map group
        mmapProcComm = worldProcsComm.split((nodeId*mmapsPerNode)+mmapIdLocalToNode, worldProcRank);
        mmapProcRank = mmapProcComm.getRank();

        isMmapLead = mmapProcRank == 0;
        mmapProcsWorldRanks = new int[mmapProcsCount];
        mmapLeadWorldRankLocalToNode =
            isMmapLead
                ? worldProcRankLocalToNode
                : (q * mmapIdLocalToNode + (mmapIdLocalToNode < r
                                                ? mmapIdLocalToNode
                                                : r));
        mmapLeadWorldRank = worldProcRank - (worldProcRankLocalToNode
                                             - mmapLeadWorldRankLocalToNode);
        // Assumes continues ranks within a node
        for (int i = 0; i < mmapProcsCount; ++i){
            mmapProcsWorldRanks[i] = mmapLeadWorldRank +i;
        }

        // Leaders talk, their color is 0
        // Followers will get a communicator of color 1,
        // but will make sure they don't talk ha ha :)
        cgProcComm = worldProcsComm.split(isMmapLead ? 0 : 1, worldProcRank);
        cgProcRank = cgProcComm.getRank();
        cgProcsCount = cgProcComm.getSize();

        parallelPattern =
            "---------------------------------------------------------\n"
            + "Machine:" + machineName + ' ' + threadCount + 'x'
            + worldProcsPerNode + 'x' + nodeCount;
        utils.printMessage(parallelPattern);
    }

    private static int[] findQandR() throws MPIException {
        int q,r;
        String str = worldProcRank+ "@" +machineName +'#';
        intBuffer.put(0, str.length());
        worldProcsComm.allReduce(intBuffer, 1, MPI.INT, MPI.MAX);
        int maxLength = intBuffer.get(0);
        CharBuffer buffer = MPI.newCharBuffer(maxLength*worldProcsCount);
        buffer.position(maxLength*worldProcRank);
        buffer.put(str);
        for (int i = str.length(); i < maxLength; ++i){
            buffer.put(i, '~');
        }

        worldProcsComm.allGather(buffer, maxLength, MPI.CHAR);
        buffer.position(0);
        Pattern nodeSep = Pattern.compile("#~*");
        Pattern nameSep = Pattern.compile("@");
        String[] nodeSplits = nodeSep.split(buffer.toString());
        HashMap<String, Integer> nodeToProcCount = new HashMap<>();
        HashMap<Integer, String> rankToNode = new HashMap<>();
        String node;
        int rank;
        String[] splits;
        for(String s: nodeSplits){
            splits = nameSep.split(s);
            rank = Integer.parseInt(splits[0].trim());
            node = splits[1].trim();
            if (nodeToProcCount.containsKey(node)){
                nodeToProcCount.put(node, nodeToProcCount.get(node)+1);
            } else {
                nodeToProcCount.put(node, 1);
            }
            rankToNode.put(rank, node);
        }

        // The following logic assumes MPI ranks are continuous within a node
        String myNode = rankToNode.get(worldProcRank);
        HashSet<String> visited = new HashSet<>();
        int rankOffset=0;
        nodeId = 0;
        for (int i = 0; i < worldProcRank; ++i){
            node = rankToNode.get(i);
            if (visited.contains(node)) continue;
            visited.add(node);
            if (node.equals(myNode)) break;
            ++nodeId;
            rankOffset += nodeToProcCount.get(node);
        }
        worldProcRankLocalToNode = worldProcRank - rankOffset;
        final int procCountOnMyNode = nodeToProcCount.get(myNode);
        q = procCountOnMyNode / mmapsPerNode;
        r = procCountOnMyNode % mmapsPerNode;

        return new int[]{q,r};
    }

    public static void tempBreak() throws MPIException {
        if (worldProcRank ==0){
            MPI.Finalize();
            System.exit(0);
        } else {
            MPI.Finalize();
            System.exit(0);
        }
    }

    public static void tearDownParallelism() throws MPIException {
        // End MPI
        MPI.Finalize();
    }

    public static void setParallelDecomposition(int globalRowCount, int targetDimension)
        throws IOException, MPIException {

        threadComm = new ThreadCommunicator(threadCount, globalRowCount, targetDimension);
        //	First divide points among processes
        procRowRanges = RangePartitioner.partition(globalRowCount,
                                                       worldProcsCount);
        Range rowRange = procRowRanges[worldProcRank]; // The range of points for this process

        procRowRange = rowRange;
        procRowStartOffset = rowRange.getStartIndex();
        procRowCount = rowRange.getLength();
        globalColCount = globalRowCount;

        // Next partition points per process among threads
        threadRowRanges = RangePartitioner.partition(procRowCount, threadCount);
        threadRowCounts = new int[threadCount];
        threadRowStartOffsets = new int[threadCount];
        IntStream.range(0, threadCount)
            .parallel()
            .forEach(threadIdx -> {
                         Range threadRowRange = threadRowRanges[threadIdx];
                         threadRowCounts[threadIdx] =
                             threadRowRange.getLength();
                         threadRowStartOffsets[threadIdx] =
                             threadRowRange.getStartIndex();
                     });

        // Allocate timing buffers
        mpiOnlyBuffer = MPI.newLongBuffer(worldProcsCount);
        threadsAndMPIBuffer = MPI.newLongBuffer(worldProcsCount * threadCount);

        cgProcsMmapRowCounts = new int[cgProcsCount];
        cgProcsMmapXByteExtents = new int[cgProcsCount];
        cgProcsMmapXDisplas = new int[cgProcsCount];

        mmapProcsRowCount = IntStream.range(mmapLeadWorldRank,
                                            mmapLeadWorldRank + mmapProcsCount)
            .map(i -> procRowRanges[i].getLength())
            .sum();
        if (isMmapLead){
            cgProcsMmapRowCounts[cgProcRank] = mmapProcsRowCount;
            cgProcComm.allGather(cgProcsMmapRowCounts, 1, MPI.INT);
            for (int i = 0; i < cgProcsCount; ++i){
                cgProcsMmapXByteExtents[i] = cgProcsMmapRowCounts[i] * targetDimension * Double.BYTES;
            }

            cgProcsMmapXDisplas[0] = 0;
            System.arraycopy(cgProcsMmapXByteExtents, 0, cgProcsMmapXDisplas, 1, cgProcsCount - 1);
            Arrays.parallelPrefix(cgProcsMmapXDisplas, (m, n) -> m + n);
        }

        boolean status = new File(mmapScratchDir).mkdirs();

        // Generate UUID and bcast it
        byte[] bytes = null;
        if (worldProcRank == 0){
            UUID id = UUID.randomUUID();
            bytes = id.toString().getBytes();
            intBuffer.put(0, bytes.length);
        }

        worldProcsComm.bcast(intBuffer, 1, MPI.INT, 0);
        int length = intBuffer.get(0);
        if (worldProcRank != 0) {
            bytes = new byte[length];
        }
        worldProcsComm.bcast(bytes, bytes.length, MPI.BYTE, 0);
        String uuid = new String(bytes);

        mmapLockFileNameOne = machineName + ".mmapId." + mmapIdLocalToNode + ".mmapLockOne." + uuid + ".bin";
        mmapEntryLockFileName = machineName + ".mmapId." + mmapIdLocalToNode + ".mmapEntryLock." + uuid + ".bin";

        final String mmapXFname = machineName + ".mmapId." + mmapIdLocalToNode + ".mmapX." + uuid +".bin";
        final String fullXFname = machineName + ".mmapId." + mmapIdLocalToNode +".fullX." + uuid + ".bin";
        try (FileChannel mmapXFc = FileChannel.open(Paths.get(mmapScratchDir,
                                                              mmapXFname),
                                                    StandardOpenOption
                                                        .CREATE,
                                                    StandardOpenOption.READ,
                                                    StandardOpenOption
                                                        .WRITE);
            FileChannel fullXFc = FileChannel.open(Paths.get(mmapScratchDir,
                                                             fullXFname),
                                                   StandardOpenOption.CREATE,StandardOpenOption.WRITE,StandardOpenOption.READ)) {


            int mmapXReadByteExtent = mmapProcsRowCount * targetDimension * Double.BYTES;
            long mmapXReadByteOffset = 0L;
            int mmapXWriteByteExtent = procRowCount * targetDimension * Double.BYTES;
            long
                mmapXWriteByteOffset =
                (procRowStartOffset - procRowRanges[mmapLeadWorldRank].getStartIndex())
                * targetDimension * Double.BYTES;
            int fullXByteExtent = globalRowCount * targetDimension * Double.BYTES;
            long fullXByteOffset = 0L;

            mmapXReadBytes = ByteBufferBytes.wrap(mmapXFc.map(
                FileChannel.MapMode.READ_WRITE, mmapXReadByteOffset,
                mmapXReadByteExtent));
            mmapXReadByteBuffer = mmapXReadBytes.sliceAsByteBuffer(
                mmapXReadByteBuffer);

            mmapXReadBytes.position(0);
            mmapXWriteBytes = mmapXReadBytes.slice(mmapXWriteByteOffset,
                                                   mmapXWriteByteExtent);

            fullXBytes = ByteBufferBytes.wrap(fullXFc.map(FileChannel.MapMode
                                                              .READ_WRITE,
                                                          fullXByteOffset,
                                                          fullXByteExtent));
            fullXByteBuffer = fullXBytes.sliceAsByteBuffer(fullXByteBuffer);

            File lockFile = new File(mmapScratchDir, mmapLockFileNameOne);
            FileChannel fc = new RandomAccessFile(lockFile, "rw").getChannel();
            mmapLockOne = ByteBufferBytes.wrap(fc.map(FileChannel.MapMode.READ_WRITE, 0, 64));
            if (isMmapLead){
                mmapLockOne.writeBoolean(FLAG, false);
                mmapLockOne.writeLong(COUNT, 0);
            }

            lockFile = new File(mmapScratchDir, mmapEntryLockFileName);
            fc = new RandomAccessFile(lockFile, "rw").getChannel();
            mmapEntryLock = ByteBufferBytes.wrap(fc.map(FileChannel.MapMode.READ_WRITE, 0, 64));
            if (isMmapLead){
                mmapEntryLock.writeBoolean(FLAG, false);
                mmapEntryLock.writeLong(COUNT, 0);
            }
        }

        /* Allocate memory maps for single double valued communications like AllReduce */
        final String mmapSFname = machineName + ".mmapId." + mmapIdLocalToNode + ".mmapS.bin";
        try (FileChannel mmapSFc = FileChannel
            .open(Paths.get(mmapScratchDir, mmapSFname),
                StandardOpenOption.CREATE, StandardOpenOption.READ,
                StandardOpenOption.WRITE)) {

            int mmapSReadByteExtent = mmapProcsCount * Double.BYTES;
            long mmapSReadByteOffset = 0L;
            int mmapSWriteByteExtent = Double.BYTES;
            long mmapSWriteByteOffset = mmapProcRank * Double.BYTES;


            mmapSReadBytes = ByteBufferBytes.wrap(mmapSFc.map(
                FileChannel.MapMode.READ_WRITE, mmapSReadByteOffset,
                mmapSReadByteExtent));
            mmapSReadByteBuffer = mmapSReadBytes.sliceAsByteBuffer(
                                                    mmapSReadByteBuffer);

            mmapSReadBytes.position(0);
            mmapSWriteBytes = mmapSReadBytes.slice(mmapSWriteByteOffset,
                mmapSWriteByteExtent);
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

    public static void allGather() throws MPIException {
        /* Safety logic to make sure all procs in the mmap has reached here. Otherwise, it's possible that
            * one (or more) procs from a same mmap may have come here while a previous call to this collective
            * is being carried out by the other procs in the same mmap. Also, note the use of a separate lock for this,
            * without that there's a chance to crash/hang */
        if (mmapEntryLock.addAndGetInt(COUNT, 1) == mmapProcsCount){
            mmapEntryLock.writeInt(COUNT, 0);
        } else {
            int count;
            do  {
                count = mmapEntryLock.readInt(COUNT);
            } while (count != 0);
        }

        if (ParallelOps.isMmapLead) {
            cgProcComm.allGatherv(mmapXReadByteBuffer,
                    cgProcsMmapXByteExtents[cgProcRank], MPI.BYTE,
                    fullXByteBuffer, cgProcsMmapXByteExtents,
                    cgProcsMmapXDisplas, MPI.BYTE);

            if (mmapProcsCount > 1) {
                mmapLockOne.writeInt(COUNT, 1); // order matters as no locks
                mmapLockOne.writeBoolean(FLAG, true);
            } else {
                /* This is for the case if you only have 1 proc per mmap,
                * then it needs to clear the flag and reset the count.
                * We special case when 1 proc per mmap under uniform mode, but
                * in a heterogeneous setting it's possible to have an mmap with 1 proc, hence this logic*/
                mmapLockOne.writeInt(COUNT, 0); // order does NOT matter for this case
                mmapLockOne.writeBoolean(FLAG, false);
            }
        } else {
            busyWaitTillDataReady();
        }
    }

    private static void busyWaitTillDataReady(){
        boolean ready = false;
        int count;
        while (!ready){
            ready = mmapLockOne.readBoolean(FLAG);
        }
        count = mmapLockOne.addAndGetInt(COUNT,1);
        if (count == mmapProcsCount){
            mmapLockOne.writeBoolean(FLAG, false);
            mmapLockOne.writeInt(COUNT, 0);
        }
    }

    public static void partialXAllGather() throws MPIException {
        cgProcComm.allGatherv(mmapXReadByteBuffer,
                              cgProcsMmapXByteExtents[cgProcRank], MPI.BYTE,
                              fullXByteBuffer, cgProcsMmapXByteExtents,
                              cgProcsMmapXDisplas, MPI.BYTE);
    }

    public static void partialSAllReduce(Op op) throws MPIException{
        cgProcComm.allReduce(mmapSReadByteBuffer, 1, MPI.DOUBLE,op);
    }

    public static void broadcast(ByteBuffer buffer, int extent, int root)
        throws MPIException {
        worldProcsComm.bcast(buffer, extent, MPI.BYTE, root);
    }

    public static void gather(LongBuffer buffer, int count, int root)
        throws MPIException {
        worldProcsComm.gather(buffer, count, MPI.LONG, root);
    }
}
