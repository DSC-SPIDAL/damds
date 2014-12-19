package cgl.imr.samples.dacidr.wdasmacof.vary;

/*
 * @author Yang Ruan(yangruan@indiana.edu)
 */
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;

public class MDSShortDataSplit {
	public static void main(String[] args) {
		if (args.length != 8) {
			System.out.println("Usage: ");
			System.out.println("[1. Data File ]");
			System.out.println("[2. Temporary directory to split data ]");
			System.out.println("[3. Temp file prefix ]");
			System.out.println("[4. Output IDs file ]");
			System.out.println("[5. Num map tasks ]");
			System.out.println("[6. row size ]");
			System.out.println("[7. column size]");
			System.out.println("[8. Type of input value format (0: short; 1: double)]");
			System.exit(0);
		}
		double beginTime = System.currentTimeMillis();

		String dataFile = args[0];
		String tmpDir = args[1];
		String tmpFilePrefix = args[2];
		String idsFile = args[3];
		int numMapTasks = Integer.valueOf(args[4]);
		int row = Integer.valueOf(args[5]);
		int col = Integer.valueOf(args[6]);
		int choice = Integer.parseInt(args[7]);

		// Create a temporary directory to hold data splits.
		if (!(new File(tmpDir)).exists()) {
			if (!(new File(tmpDir)).mkdir()) {
				System.err
						.print("Failed to create the temporary directory to split data");
				System.exit(-1);
			}
		}
		try {
			if(choice == 0)
				splitDataforShort(dataFile, tmpDir, tmpFilePrefix, idsFile, numMapTasks, row, col);
			else if (choice == 1)
				splitDataforDouble(dataFile, tmpDir, tmpFilePrefix, idsFile, numMapTasks, row, col);
			else{
				System.err.println("The choice must be 1 or 0");
				System.exit(2);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		double endTime = System.currentTimeMillis();
		System.out
				.println("==========================================================");
		System.out.println("Time to split data = " + (endTime - beginTime)
				/ 1000 + " Seconds.");
		System.out
				.println("==========================================================");
		System.exit(0);
	}

	private static void splitDataforDouble(String dataFile, String tmpDir,
			String tmpFilePrefix, String idsFile, int numMapTasks, int row, int col)
			throws IOException {
		BufferedInputStream reader = new BufferedInputStream(new FileInputStream(
				dataFile));
		DataInputStream din = new DataInputStream(reader);
		BufferedWriter idsWriter = new BufferedWriter(new FileWriter(idsFile));

		String outputFile = null;
		int blockHeight = row / numMapTasks;
		int rem = row % numMapTasks;

		double[][] rows;
		int start = 0;
		int end = 0;
		int curHeight = 0;
		int count = 0;

		for (int i = 0; i < numMapTasks; i++) {
			System.out.println("The " + i + "th maptask");
			idsWriter.write(i + "\t");
			outputFile = tmpDir + "/" + tmpFilePrefix + i;

			end += blockHeight;
			if (rem > 0) {
				end++;
				rem--;
			}
			curHeight = end - start;
			rows = new double[curHeight][col];
			count = 0;
			for (int j = start; j < end; j++) {
				for (int k = 0; k < col; k++) { // Note the start
					rows[count][k] = din.readDouble();
				}
				count++;
			}
			writeToBinFile(rows, curHeight, col, outputFile);
			idsWriter.write(curHeight + "\t" + col + "\t" + i + "\t" + start + "\n");
			start = end;
		}
		din.close();
		idsWriter.close();
	}
	
	public static void writeToBinFile(double[][] row, int curHeight, int size, String fileName) throws IOException {
		BufferedOutputStream bout = new BufferedOutputStream(
				new FileOutputStream(fileName));
		DataOutputStream dout = new DataOutputStream(bout);

		// First two parameters are the dimensions.
		//dout.writeInt(height);
		//dout.writeInt(width);
		//dout.writeInt(row);
		//dout.writeInt(rowOffset);
		for (int i = 0; i < curHeight; i++) {
			for (int j = 0; j < size; j++) {
				dout.writeDouble(row[i][j]);
			}
		}
		dout.flush();
		bout.flush();
		dout.close();
		bout.close();
	}
	
	/**
	 * data size should be provided for reading from bin file
	 * 
	 * @param dataFile
	 * @param tmpDir
	 * @param tmpFilePrefix
	 * @param idsFile
	 * @param numMapTasks
	 * @param size
	 * @throws IOException
	 */
	
	private static void splitDataforShort(String dataFile, String tmpDir,
			String tmpFilePrefix, String idsFile, int numMapTasks, int row, int col)
			throws IOException {
		BufferedInputStream reader = new BufferedInputStream(new FileInputStream(
				dataFile));
		DataInputStream din = new DataInputStream(reader);
		BufferedWriter idsWriter = new BufferedWriter(new FileWriter(idsFile));

		String outputFile = null;
		int blockHeight = row / numMapTasks;
		int rem = row % numMapTasks;

		MDSShortMatrixData rowData;
		short[][] rows;
		int start = 0;
		int end = 0;
		int curHeight = 0;
		int count = 0;

		for (int i = 0; i < numMapTasks; i++) {
			System.out.println("The " + i + "th maptask");
			idsWriter.write(i + "\t");
			outputFile = tmpDir + "/" + tmpFilePrefix + i;

			end += blockHeight;
			if (rem > 0) {
				end++;
				rem--;
			}
			curHeight = end - start;
			rows = new short[curHeight][col];
			count = 0;
			for (int j = start; j < end; j++) {
				for (int k = 0; k < col; k++) { // Note the start
					rows[count][k] = din.readShort();
				}
				count++;
			}
			rowData = new MDSShortMatrixData(rows, curHeight, col, i, start);
			idsWriter.write(curHeight + "\t" + col + "\t" + i + "\t" + start + "\n");
			start = end;
			rowData.writeToBinFile(outputFile);
		}
		din.close();
		idsWriter.close();
	}
}
