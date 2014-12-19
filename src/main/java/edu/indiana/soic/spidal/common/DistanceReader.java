package edu.indiana.soic.spidal.common;

import com.google.common.io.LittleEndianDataInputStream;

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

public interface DistanceReader {
    public double getDistance(int globalRow, int globalCol);

    public static DistanceReader readRowRange(String fname, Range rows, int globalColCount, ByteOrder
            endianness, boolean mmap){
        if (mmap) {
            try (FileChannel fc = (FileChannel) Files.newByteChannel(Paths.get(fname), StandardOpenOption.READ)) {
                int dataTypeSize = Short.BYTES;
                long pos = ((long) rows.getStartIndex()) * globalColCount * dataTypeSize;
                MappedByteBuffer mappedBytes = fc.map(FileChannel.MapMode.READ_ONLY, pos,
                        rows.getLength() * globalColCount * dataTypeSize);
                mappedBytes.order(endianness);
                return (globalRow, globalCol) -> {
                    int position = (globalRow - rows.getStartIndex()) * globalColCount + globalCol; // element position - not the byte position
                    return mappedBytes.getShort(position * dataTypeSize) / (Short.MAX_VALUE * 1.0); // pos*2 is the byte position
                };
            } catch (IOException e) {
                e.printStackTrace();
            }
        } else {
            int startRow = rows.getStartIndex();
            int numRows = rows.getLength();
            try (FileInputStream fis = new FileInputStream(fname)) {
                DataInput di = endianness == ByteOrder.BIG_ENDIAN ? new DataInputStream(
                        fis) : new LittleEndianDataInputStream(fis);

                int numBytesToSkip = startRow * globalColCount * Short.BYTES;
                int skippedBytes = di.skipBytes(numBytesToSkip);
                if (skippedBytes != numBytesToSkip)
                    throw new IOException(String.format("Requested %1$d bytes to skip, but could skip only %2$d bytes", numBytesToSkip, skippedBytes));

                short[][] buffer = new short[numRows][];
                for (int i = 0; i < numRows; ++i) {
                    buffer[i] = new short[globalColCount];
                    for (int j = 0; j < globalColCount; ++j) {
                        buffer[i][j] = di.readShort();
                    }
                }
                return (globalRow, globalCol) -> {
                    int localRow = globalRow - rows.getStartIndex();
                    return buffer[localRow][globalCol] / (1.0*Short.MAX_VALUE);
                };
            } catch (IOException e){
                e.printStackTrace();
            }
        }
        return null;
    }

}

