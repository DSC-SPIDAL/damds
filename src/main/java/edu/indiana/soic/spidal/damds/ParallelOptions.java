package edu.indiana.soic.spidal.damds;

import edu.indiana.soic.spidal.common.Range;
import edu.indiana.soic.spidal.common.RangePartitioner;
import mpi.Intracomm;
import mpi.MPI;
import mpi.MPIException;

import static edu.rice.hj.Module0.finalizeHabanero;
import static edu.rice.hj.Module0.initializeHabanero;

public class ParallelOptions {
    public static int nodeCount=1;
    public static int threadCount=1;

    public static Intracomm mpiComm;
    public static int mpiRank;
    public static int mpiSize;
    public static String parallelPattern;

    public static int localRowStartOffset;
    public static int localRowCount;
    public static Range localRowRange;

    public static void setupParallelism(String[] args) throws MPIException {
        // Set up threads
        initializeHabanero();

        //  Set up MPI
        MPI.Init(args);
        mpiComm = MPI.COMM_WORLD; //initializing MPI world communicator
        mpiRank = mpiComm.getRank();
        mpiSize = mpiComm.getSize();

        // Set up MPI
        int mpiPerNode = mpiSize / nodeCount;

        if ((mpiPerNode * nodeCount) != mpiSize)
        {
            Utils.printAndThrowRuntimeException("Inconsistent MPI counts Nodes " + nodeCount + " Size " + mpiSize);
        }

        parallelPattern = "---------------------------------------------------------\nMachine:" + MPI.getProcessorName() + " " + threadCount + "x" + mpiPerNode + "x" + nodeCount;
        Utils.printMessage(parallelPattern);
    }

    public static void tearDownParallelism() throws MPIException {
        // Finalize threads
        finalizeHabanero();

        // End MPI
        MPI.Finalize();
    }

    public static void setParallelDecomposition(int globalRowCount) {
        //	First divide points among processes
        Range[] rowRanges = RangePartitioner.Partition(globalRowCount, mpiSize);
        Range rowRange = rowRanges[mpiRank]; // The range of points for this process

        localRowRange = rowRange;
        localRowStartOffset = rowRange.getStartIndex();
        localRowCount = rowRange.getLength();
    }
}
