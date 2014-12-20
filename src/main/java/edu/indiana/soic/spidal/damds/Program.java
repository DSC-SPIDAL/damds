package edu.indiana.soic.spidal.damds;

import com.google.common.base.Optional;
import edu.indiana.soic.spidal.common.BinaryReader;
import edu.indiana.soic.spidal.common.DoubleStatistics;
import edu.indiana.soic.spidal.configuration.ConfigurationMgr;
import edu.indiana.soic.spidal.configuration.section.DAMDSSection;
import mpi.MPIException;
import org.apache.commons.cli.*;

import java.nio.ByteOrder;
import java.util.DoubleSummaryStatistics;
import java.util.stream.IntStream;

public class Program {
    private static Options programOptions = new Options();
    static {
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_C),Constants.CMD_OPTION_LONG_C, true,
                Constants.CMD_OPTION_DESCRIPTION_C);
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_N),Constants.CMD_OPTION_LONG_N,true,
                Constants.CMD_OPTION_DESCRIPTION_N);
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_T),Constants.CMD_OPTION_LONG_T,true,
                Constants.CMD_OPTION_DESCRIPTION_T);
    }

    //Config Settings
    public static DAMDSSection config;
    public static ByteOrder byteOrder;
    public static BinaryReader distances;
    public static BinaryReader weights;


    /**
     * Weighted SMACOF based on Deterministic Annealing algorithm
     * @param args command line arguments to the program, which should include
     *             -c "path to config file" -t "number of threads" -n "number of nodes"
     *             The options may also be given as longer names
     *             --configFile, --threadCount, and --nodeCount respectively
     */
    public static void main(String[] args) {
        Optional<CommandLine> parserResult = parseCommandLineArguments(args, programOptions);
        if (!parserResult.isPresent()){
            System.out.println(Constants.ERR_PROGRAM_ARGUMENTS_PARSING_FAILED);
            new HelpFormatter().printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        CommandLine cmd = parserResult.get();
        if (!(cmd.hasOption(Constants.CMD_OPTION_LONG_C) &&
                cmd.hasOption(Constants.CMD_OPTION_LONG_N) &&
                cmd.hasOption(Constants.CMD_OPTION_LONG_T))){
            System.out.println(Constants.ERR_INVALID_PROGRAM_ARGUMENTS);
            new HelpFormatter().printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        //  Read Metadata using this as source of other metadata
        ReadControlFile(cmd);

        try {
            //  Set up MPI and threads parallelism
            ParallelOps.setupParallelism(args);
            ParallelOps.setParallelDecomposition(config.numberDataPoints);
            distances = BinaryReader
                    .readRowRange(config.distanceMatrixFile, ParallelOps.localRowRange, ParallelOps.globalColCount,
                                  byteOrder, config.isMemoryMapped, true);
            DoubleStatistics distanceSummary;
            if (!config.isSammon) {
                weights = BinaryReader
                        .readRowRange(config.weightMatrixFile, ParallelOps.localRowRange, ParallelOps.globalColCount,
                                      byteOrder, config.isMemoryMapped, false);
                // Non Sammon mode - use only distances that have non zero corresponding weights
                distanceSummary =
                        IntStream.range(0, ParallelOps.localRowCount * ParallelOps.globalColCount).parallel()
                                 .filter(i -> weights
                                         .getValue(i / ParallelOps.globalColCount, i % ParallelOps.globalColCount) != 0)
                                 .mapToDouble(i -> distances.getValue(i / ParallelOps.globalColCount,
                                                                      i % ParallelOps.globalColCount))
                                 .collect(DoubleStatistics::new, DoubleStatistics::accept, DoubleStatistics::combine);
            } else {
                // Sammon mode - use all distances
                distanceSummary =
                        IntStream.range(0, ParallelOps.localRowCount * ParallelOps.globalColCount).parallel()
                                 .mapToDouble(i -> distances
                                         .getValue(i / ParallelOps.globalColCount, i % ParallelOps.globalColCount))
                                 .collect(DoubleStatistics::new, DoubleStatistics::accept, DoubleStatistics::combine);
            }
            distanceSummary = ParallelOps.allReduce(distanceSummary);
            Utils.printMessage(distanceSummary.toString());


        } catch (MPIException e) {
            Utils.printAndThrowRuntimeException(new RuntimeException(e));
        }



    }

    private static void ReadControlFile(CommandLine cmd) {
        config = ConfigurationMgr.LoadConfiguration(cmd.getOptionValue(Constants.CMD_OPTION_LONG_C)).damdsSection;
        ParallelOps.nodeCount = Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_N));
        ParallelOps.threadCount = Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_T));
        byteOrder = config.isBigEndian ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN;
    }

    /**
     * Parse command line arguments
     * @param args Command line arguments
     * @param opts Command line options
     * @return An <code>Optional&lt;CommandLine&gt;</code> object
     */
    private static Optional<CommandLine> parseCommandLineArguments(String [] args, Options opts){

        CommandLineParser optParser = new GnuParser();

        try {
            return Optional.fromNullable(optParser.parse(opts, args));
        } catch (ParseException e) {
            System.out.println(e);
        }
        return Optional.fromNullable(null);
    }
}
