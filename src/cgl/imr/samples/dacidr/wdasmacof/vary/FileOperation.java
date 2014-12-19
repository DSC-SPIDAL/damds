package cgl.imr.samples.dacidr.wdasmacof.vary;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;

public class FileOperation {
	public static int[][] loadV(String fileName, int row, int col) throws IOException{
		DataInputStream din = new DataInputStream(new BufferedInputStream(new FileInputStream(fileName)));
		int[][] VPlus = new int[row][col];
		for(int i = 0; i < row; ++i)
			for(int j = 0; j < col; ++j)
				VPlus[i][j] = (int) din.readDouble();
		din.close();
		return VPlus;
	}
	
	public static double[][] loadVPlus(String fileName, int row, int col) throws IOException{
		DataInputStream din = new DataInputStream(new BufferedInputStream(new FileInputStream(fileName)));
		double[][] VPlus = new double[row][col];
		for(int i = 0; i < row; ++i)
			for(int j = 0; j < col; ++j)
				VPlus[i][j] = din.readDouble();
		din.close();
		return VPlus;
	}
	
	public static short[][] loadWeights(String fileName, int row, int col) throws IOException{
		DataInputStream din = new DataInputStream(new BufferedInputStream(new FileInputStream(fileName)));
		short[][] weights = new short[row][col];
		for(int i = 0; i < row; ++i)
			for(int j = 0; j < col; ++j)
				weights[i][j] = din.readShort();
		din.close();
		return weights;
	}

    /**
     * Creates weights for Sammon mapping
     * @param originalDistances distances in the original space as shorts.
     *                          These are expected to in the range [0,1] and are represented as shorts by
     *                          multiplying by Short.MAX_VALUE
     * @param originalMeanDistance mean (average) distance in the original space out of all distances.
     *                             Expected to be in the range [0,1]
     * @param row number of rows
     * @param col number of cols
     * @return a two dimensional short array consisting weights. The formula for weights is
     * W_ij = 1/Max(d_ij, 0.001*meanD)
     * where _ij means index, W means weight, d means original distance, and meanD is originalMeanDistance
     */
    public static short[][] loadSammonWeights(short [][] originalDistances, double originalMeanDistance, int row, int col){
        short[][] weights = new short[row][col];
        for(int i = 0; i < row; ++i)
            for(int j = 0; j < col; ++j)
                weights[i][j] = (short)((1.0/Math.max(originalDistances[i][j]*1.0/Short.MAX_VALUE, 0.001 * originalMeanDistance)) * Short.MAX_VALUE);
        return weights;
    }
}
