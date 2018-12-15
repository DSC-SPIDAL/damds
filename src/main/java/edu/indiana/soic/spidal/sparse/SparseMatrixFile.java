package edu.indiana.soic.spidal.sparse;

import edu.indiana.soic.spidal.common.Range;
import org.apache.commons.lang3.ArrayUtils;

import java.io.*;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SparseMatrixFile {



    private String pathname;

    private File indicesFile;
    private File dataFile;



    private SparseMatrixFile() {
    }

    public SparseMatrixFile(String pathname, long rowLength, long colLength) {
        this.pathname = pathname;
    }

    /**
     * Square matrix only
     *
     * @param indicesPath indices
     * @param dataPath    data
     * @return partial sparse matrix
     */
    public static SparseMatrix loadIntoMemory(String indicesPath, String dataPath, Range globalThreadRowRange, int dim) {
        //TODO: check if we can use arrays instead of lists we need to know the length of values before hand for this
        // maybe we can do a two pass method
        int startRow = globalThreadRowRange.getStartIndex();
        int endRow = globalThreadRowRange.getEndIndex();
        int length = globalThreadRowRange.getLength();

        if (startRow < 0 || startRow > endRow || startRow > dim) {
            throw new RuntimeException("Illegal row range");
        }
        try {
            SparseMatrixFile smf = new SparseMatrixFile();
            smf.indicesFile = new File(indicesPath);
            smf.dataFile = new File(dataPath);
            BufferedInputStream indexIn = new BufferedInputStream(new FileInputStream(smf.indicesFile));
            BufferedInputStream dataIn = new BufferedInputStream(new FileInputStream(smf.dataFile));
            byte[] buf = new byte[8];
            int len = 0;
            List<Double> values = new ArrayList<>();
            List<Integer> columns = new ArrayList<>();
            int[] rowPointer = new int[length];
            Arrays.fill(rowPointer, -1);

            int count = 0;
            while (indexIn.available() > 0 && dataIn.available() > 0) {
                len = indexIn.read(buf);
                long i = bytesToLong(buf);
                len = indexIn.read(buf);
                long j = bytesToLong(buf);
                len = dataIn.read(buf);
                double value = ByteBuffer.wrap(buf).getDouble();
                if (i >= startRow && i <= endRow) {
                    int localRow = (int) i - startRow;
                    values.add(value);
                    columns.add((int)j);
                    if(rowPointer[localRow] == -1){
                        rowPointer[localRow] = count;
                    }
                    count++;
                }
            }
            SparseMatrix sparseMatrix = new SparseMatrix(ArrayUtils.toPrimitive((Double[])values.toArray()),
                    ArrayUtils.toPrimitive((Integer[])columns.toArray()), rowPointer);

            indexIn.close();
            dataIn.close();
            return sparseMatrix;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    private static long bytesToLong(byte[] bytes) {
        return (bytes[0] & 0xFFL) << 56
                | (bytes[1] & 0xFFL) << 48
                | (bytes[2] & 0xFFL) << 40
                | (bytes[3] & 0xFFL) << 32
                | (bytes[4] & 0xFFL) << 24
                | (bytes[5] & 0xFFL) << 16
                | (bytes[6] & 0xFFL) << 8
                | (bytes[7] & 0xFFL);
    }

    private static byte[] longToBytes(long l) {
        byte[] result = new byte[8];
        for (int i = 7; i >= 0; i--) {
            result[i] = (byte) (l & 0xFF);
            l >>= 8;
        }
        return result;
    }
}
