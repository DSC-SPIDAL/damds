package edu.indiana.soic.spidal.damds;

public class Constants {
    static final String PROGRAM_NAME = "DAMDS";

    static final char CMD_OPTION_SHORT_C = 'c';
    static final String CMD_OPTION_LONG_C = "configFile";
    static final String CMD_OPTION_DESCRIPTION_C = "Configuration file";
    static final char CMD_OPTION_SHORT_N = 'n';
    static final String CMD_OPTION_LONG_N = "nodeCount";
    static final String CMD_OPTION_DESCRIPTION_N = "Node count";
    static final char CMD_OPTION_SHORT_T = 't';
    static final String CMD_OPTION_LONG_T = "threadCount";
    public static final String CMD_OPTION_SHORT_MMAPS = "mmaps";
    public static final String
        CMD_OPTION_DESCRIPTION_MMAPS =
        "Number of memory mapped groups per node";
    public static final String CMD_OPTION_SHORT_MMAP_SCRATCH_DIR = "mmapdir";
    public static final String
        CMD_OPTION_DESCRIPTION_MMAP_SCRATCH_DIR =
        "Scratch directory to store memmory mapped files. A node local "
        + "volatile storage like tmpfs is advised for this";

    static final String CMD_OPTION_DESCRIPTION_T = "Thread count";
    static final String
        ERR_PROGRAM_ARGUMENTS_PARSING_FAILED =
        "Argument parsing failed!";
    static final String
        ERR_INVALID_PROGRAM_ARGUMENTS =
        "Invalid program arguments!";
    static final String ERR_EMPTY_FILE_NAME = "File name is null or empty!";

    public static String errWrongNumOfBytesSkipped(
        int requestedBytesToSkip, int numSkippedBytes) {
        String
            msg =
            "Requested %1$d bytes to skip, but could skip only %2$d bytes";
        return String.format(msg, requestedBytesToSkip, numSkippedBytes);
    }
}
