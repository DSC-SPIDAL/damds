package edu.indiana.soic.spidal.damds;

import com.google.common.base.Optional;
import edu.indiana.soic.spidal.common.DistanceReader;
import edu.indiana.soic.spidal.configuration.ConfigurationMgr;
import edu.indiana.soic.spidal.configuration.section.DAMDSSection;
import mpi.MPIException;
import org.apache.commons.cli.*;

import java.nio.ByteOrder;

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
    public static DistanceReader distances;


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
            ParallelOptions.setupParallelism(args);
            ParallelOptions.setParallelDecomposition(config.numberDataPoints);
            distances = DistanceReader.readRowRange(config.distanceMatrixFile,ParallelOptions.localRowRange, config.numberDataPoints, Short.BYTES, byteOrder, config.isMemoryMapped);

        } catch (MPIException e) {
            Utils.printAndThrowRuntimeException(new RuntimeException(e));
        }



    }

    private static void ReadControlFile(CommandLine cmd) {
        config = ConfigurationMgr.LoadConfiguration(cmd.getOptionValue(Constants.CMD_OPTION_LONG_C)).damdsSection;
        ParallelOptions.nodeCount = Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_N));
        ParallelOptions.threadCount = Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_T));
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
