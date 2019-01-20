package edu.indiana.soic.spidal.damds;

import com.google.common.base.Optional;
import com.google.common.base.Stopwatch;
import edu.indiana.soic.spidal.configuration.ConfigurationMgr;
import edu.indiana.soic.spidal.configuration.section.DAMDSSection;
import edu.indiana.soic.spidal.damds.threads.SpidalThreads;
import edu.indiana.soic.spidal.damds.threads.ThreadBitAssigner;
import mpi.MPIException;
import net.openhft.affinity.Affinity;
import org.apache.commons.cli.*;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.ByteOrder;
import java.nio.LongBuffer;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Date;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;


public class SparseProgram {
    private static Options programOptions = new Options();

    static {
        programOptions.addOption(
            String.valueOf(Constants.CMD_OPTION_SHORT_C),
            Constants.CMD_OPTION_LONG_C, true,
            Constants.CMD_OPTION_DESCRIPTION_C);
        programOptions.addOption(
            String.valueOf(Constants.CMD_OPTION_SHORT_N),
            Constants.CMD_OPTION_LONG_N, true,
            Constants.CMD_OPTION_DESCRIPTION_N);
        programOptions.addOption(
            String.valueOf(Constants.CMD_OPTION_SHORT_T),
            Constants.CMD_OPTION_LONG_T, true,
            Constants.CMD_OPTION_DESCRIPTION_T);

        programOptions.addOption(Constants.CMD_OPTION_SHORT_MMAPS, true, Constants.CMD_OPTION_DESCRIPTION_MMAPS);
        programOptions.addOption(
            Constants.CMD_OPTION_SHORT_MMAP_SCRATCH_DIR, true,
            Constants.CMD_OPTION_DESCRIPTION_MMAP_SCRATCH_DIR);
        programOptions.addOption(
                Constants.CMD_OPTION_SHORT_BIND_THREADS, true,
                Constants.CMD_OPTION_DESCRIPTION_BIND_THREADS);

        programOptions.addOption(
                Constants.CMD_OPTION_SHORT_CPS, true,
                Constants.CMD_OPTION_DESCRIPTION_CPS);


    }

    //Config Settings
    public static DAMDSSection config;
    public static ByteOrder byteOrder;

    public static int BlockSize;
    private static Utils utils = new Utils(0);
    private static boolean bind;
    private static int cps;

    private static SpidalThreads threads = null;

    /**
     * Weighted SMACOF based on Deterministic Annealing algorithm
     *
     * @param args command line arguments to the program, which should include
     *             -c path to config file
     *             -t number of threads
     *             -n number of nodes
     *             The options may also be given as longer names
     *             --configFile, --threadCount, and --nodeCount respectively
     */
    public static void  main(String[] args) {
        Optional<CommandLine> parserResult =
            parseCommandLineArguments(args, programOptions);

        if (!parserResult.isPresent()) {
            System.out.println(Constants.ERR_PROGRAM_ARGUMENTS_PARSING_FAILED);
            new HelpFormatter()
                .printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        CommandLine cmd = parserResult.get();
        if (!(cmd.hasOption(Constants.CMD_OPTION_LONG_C) &&
              cmd.hasOption(Constants.CMD_OPTION_LONG_N) &&
              cmd.hasOption(Constants.CMD_OPTION_LONG_T))) {
            System.out.println(Constants.ERR_INVALID_PROGRAM_ARGUMENTS);
            new HelpFormatter()
                .printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        try {
            //  Read Metadata using this as source of other metadata
            readConfiguration(cmd);

            //  Set up MPI and threads parallelism
            ParallelOps.setupParallelism(args);
            ParallelOps.setParallelDecomposition(
                config.numberDataPoints, config.targetDimension);

            // Note - a barrier to get cleaner timings
            ParallelOps.worldProcsComm.barrier();
            Stopwatch mainTimer = Stopwatch.createStarted();

            utils.printMessage("\n== DAMDS run started on " + new Date() + " ==\n");
            utils.printMessage(config.toString(false));

            /* TODO - Fork - join starts here */
            if (ParallelOps.threadCount > 1) {
                Lock lock = new ReentrantLock();
                launchHabaneroApp(
                    () -> forallChunked(
                        0, ParallelOps.threadCount - 1,
                        (threadIdx) -> {

                            if (bind) {
                                BitSet bitSet = ThreadBitAssigner.getBitSet(ParallelOps.worldProcRank, threadIdx, ParallelOps.threadCount, cps);
                                Affinity.setAffinity(bitSet);
                            }

                            final SparseProgramWorker worker = new SparseProgramWorker
                                    (threadIdx,
                                    ParallelOps
                                    .threadComm, config, byteOrder,
                                            BlockSize, mainTimer,
                                            lock);
                            try {
                                worker.run();
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }));
            }
            else {
                if (bind) {
                    BitSet bitSet = ThreadBitAssigner.getBitSet(ParallelOps.worldProcRank, 0, ParallelOps.threadCount, cps);
                    Affinity.setAffinity(bitSet);
                }
                new SparseProgramWorker(0, ParallelOps.threadComm, config,
                        byteOrder, BlockSize, mainTimer, null)
                        .run();
            }


            /* TODO - Fork-join should end here */


            /*
             // TODO Fix error handling here
            printTimings(totalTime, temperatureLoopTime);*/

            utils.printMessage("== DAMDS run completed on " + new Date() + " ==");

            ParallelOps.tearDownParallelism();
        }
        catch (MPIException | IOException e) {
            utils.printAndThrowRuntimeException(new RuntimeException(e));
        }
    }

    private static void readConfiguration(CommandLine cmd) {
        config = ConfigurationMgr.LoadConfiguration(
            cmd.getOptionValue(Constants.CMD_OPTION_LONG_C)).damdsSection;
        ParallelOps.nodeCount =
            Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_N));
        ParallelOps.threadCount =
            Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_T));
        ParallelOps.mmapsPerNode = cmd.hasOption(Constants.CMD_OPTION_SHORT_MMAPS) ? Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_SHORT_MMAPS)) : 1;
        ParallelOps.mmapScratchDir = cmd.hasOption(Constants.CMD_OPTION_SHORT_MMAP_SCRATCH_DIR) ? cmd.getOptionValue(Constants.CMD_OPTION_SHORT_MMAP_SCRATCH_DIR) : ".";

        byteOrder =
            config.isBigEndian ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN;
        BlockSize = config.blockSize;

        bind = !cmd.hasOption(Constants.CMD_OPTION_SHORT_BIND_THREADS) ||
                Boolean.parseBoolean(cmd.getOptionValue(Constants.CMD_OPTION_SHORT_BIND_THREADS));
        cps = (cmd.hasOption(Constants.CMD_OPTION_SHORT_CPS)) ? Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_SHORT_CPS)) : -1;
        if (cps == -1){
            utils.printMessage("Disabling thread binding as cps is not specified");
            bind = false;
        }
    }

    /**
     * Parse command line arguments
     *
     * @param args Command line arguments
     * @param opts Command line options
     * @return An <code>Optional&lt;CommandLine&gt;</code> object
     */
    private static Optional<CommandLine> parseCommandLineArguments(
        String[] args, Options opts) {

        CommandLineParser optParser = new GnuParser();

        try {
            return Optional.fromNullable(optParser.parse(opts, args));
        }
        catch (ParseException e) {
            e.printStackTrace();
        }
        return Optional.fromNullable(null);
    }
}
