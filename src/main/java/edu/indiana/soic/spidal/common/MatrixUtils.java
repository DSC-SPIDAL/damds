package edu.indiana.soic.spidal.common;

import java.util.function.ToDoubleBiFunction;

/**
 * @author Yang Ruan, yangruan@cs.indiana.edu Some code segments contain
 *         in this class are directly inherited from the C# version of the
 *         shared memory MDS program wirtten by my collegue Seung-Hee Bea.
 *         <p/>
 *         A collection of matrix related utility operations.
 */
public class MatrixUtils {
    /**
     * Block matrix multiplication.
     *
     * @param A         - Matrix
     * @param Adiag     - Diag value of A
     * @param B         - Matrix
     * @param aHeight   - height of matrix A
     * @param bWidth    - width of matrix B
     * @param comm      - size of the common face
     * @param rowOffset - offset from row index to col index
     * @return
     */
    public static double[][] matrixMultiply(
        ToDoubleBiFunction<Integer, Integer> A, double[] Adiag, double[][] B,
        int aHeight, int bWidth, int comm, int bz, int rowOffset) {
        double[][] C = new double[aHeight][bWidth];

        int aHeightBlocks = aHeight / bz; // size = Height of A
        int aLastBlockHeight = aHeight - (aHeightBlocks * bz);
        if (aLastBlockHeight > 0) {
            aHeightBlocks++;
        }

        int bWidthBlocks = bWidth / bz; // size = Width of B
        int bLastBlockWidth = bWidth - (bWidthBlocks * bz);
        if (bLastBlockWidth > 0) {
            bWidthBlocks++;
        }

        int commnBlocks = comm / bz; // size = Width of A or Height of B
        int commLastBlockWidth = comm - (commnBlocks * bz);
        if (commLastBlockWidth > 0) {
            commnBlocks++;
        }

        int aBlockHeight = bz;
        int bBlockWidth = bz;
        int commBlockWidth = bz;

        for (int ib = 0; ib < aHeightBlocks; ib++) {
            if (aLastBlockHeight > 0 && ib == (aHeightBlocks - 1)) {
                aBlockHeight = aLastBlockHeight;
            }
            bBlockWidth = bz;
            commBlockWidth = bz;
            for (int jb = 0; jb < bWidthBlocks; jb++) {
                if (bLastBlockWidth > 0 && jb == (bWidthBlocks - 1)) {
                    bBlockWidth = bLastBlockWidth;
                }
                commBlockWidth = bz;
                for (int kb = 0; kb < commnBlocks; kb++) {
                    if (commLastBlockWidth > 0 && kb == (commnBlocks - 1)) {
                        commBlockWidth = commLastBlockWidth;
                    }

                    for (int i = ib * bz; i < (ib * bz) + aBlockHeight; i++) {
                        for (int j = jb * bz; j < (jb * bz) + bBlockWidth;
                             j++) {
                            for (int k = kb * bz;
                                 k < (kb * bz) + commBlockWidth; k++) {
                                double aVal = 0;
                                if (i + rowOffset == k) {
                                    aVal = Adiag[i];
                                }
                                else {
                                    //reverse the value from weight
                                    aVal = -(A.applyAsDouble(i, k));
                                }

                                if (aVal != 0 && B[k][j] != 0) {
                                    C[i][j] += aVal * B[k][j];
                                }
                            }
                        }
                    }
                }
            }
        }

        return C;
    }

    /**
     * Block matrix multiplication.
     *
     * @param A         - Matrix
     * @param Adiag     - Diag value of A
     * @param B         - Matrix
     * @param aHeight   - height of matrix A
     * @param bWidth    - width of matrix B
     * @param comm      - size of the common face
     * @param rowOffset - offset from row index to col index
     * @return
     */
    public static double[][] matrixMultiply(
        short[][] A, double[] Adiag, double[][] B, int aHeight, int bWidth,
        int comm, int bz, int rowOffset) {
        double[][] C = new double[aHeight][bWidth];

        int aHeightBlocks = aHeight / bz; // size = Height of A
        int aLastBlockHeight = aHeight - (aHeightBlocks * bz);
        if (aLastBlockHeight > 0) {
            aHeightBlocks++;
        }

        int bWidthBlocks = bWidth / bz; // size = Width of B
        int bLastBlockWidth = bWidth - (bWidthBlocks * bz);
        if (bLastBlockWidth > 0) {
            bWidthBlocks++;
        }

        int commnBlocks = comm / bz; // size = Width of A or Height of B
        int commLastBlockWidth = comm - (commnBlocks * bz);
        if (commLastBlockWidth > 0) {
            commnBlocks++;
        }

        int aBlockHeight = bz;
        int bBlockWidth = bz;
        int commBlockWidth = bz;

        for (int ib = 0; ib < aHeightBlocks; ib++) {
            if (aLastBlockHeight > 0 && ib == (aHeightBlocks - 1)) {
                aBlockHeight = aLastBlockHeight;
            }
            bBlockWidth = bz;
            commBlockWidth = bz;
            for (int jb = 0; jb < bWidthBlocks; jb++) {
                if (bLastBlockWidth > 0 && jb == (bWidthBlocks - 1)) {
                    bBlockWidth = bLastBlockWidth;
                }
                commBlockWidth = bz;
                for (int kb = 0; kb < commnBlocks; kb++) {
                    if (commLastBlockWidth > 0 && kb == (commnBlocks - 1)) {
                        commBlockWidth = commLastBlockWidth;
                    }

                    for (int i = ib * bz; i < (ib * bz) + aBlockHeight; i++) {
                        for (int j = jb * bz; j < (jb * bz) + bBlockWidth;
                             j++) {
                            for (int k = kb * bz;
                                 k < (kb * bz) + commBlockWidth; k++) {
                                double aVal = 0;
                                if (i + rowOffset == k) {
                                    aVal = Adiag[i];
                                }
                                else {
                                    //reverse the value from weight
                                    aVal = -(A[i][k]);
                                }

                                if (aVal != 0 && B[k][j] != 0) {
                                    C[i][j] += aVal * B[k][j];
                                }
                            }
                        }
                    }
                }
            }
        }

        return C;
    }

    /**
     * Block matrix multiplication.
     *
     * @param A       - Matrix
     * @param B       - Matrix
     * @param aHeight - height of matrix A
     * @param bWidth  - width of matrix B
     * @param comm    - size of the common face
     * @return
     */
    public static double[][] matrixMultiply(
        float[][] A, double[][] B, int aHeight, int bWidth, int comm, int bz) {

        double[][] C = new double[aHeight][bWidth];

        int aHeightBlocks = aHeight / bz; // size = Height of A
        int aLastBlockHeight = aHeight - (aHeightBlocks * bz);
        if (aLastBlockHeight > 0) {
            aHeightBlocks++;
        }

        int bWidthBlocks = bWidth / bz; // size = Width of B
        int bLastBlockWidth = bWidth - (bWidthBlocks * bz);
        if (bLastBlockWidth > 0) {
            bWidthBlocks++;
        }

        int commnBlocks = comm / bz; // size = Width of A or Height of B
        int commLastBlockWidth = comm - (commnBlocks * bz);
        if (commLastBlockWidth > 0) {
            commnBlocks++;
        }

        int aBlockHeight = bz;
        int bBlockWidth = bz;
        int commBlockWidth = bz;

        for (int ib = 0; ib < aHeightBlocks; ib++) {
            if (aLastBlockHeight > 0 && ib == (aHeightBlocks - 1)) {
                aBlockHeight = aLastBlockHeight;
            }
            bBlockWidth = bz;
            commBlockWidth = bz;
            for (int jb = 0; jb < bWidthBlocks; jb++) {
                if (bLastBlockWidth > 0 && jb == (bWidthBlocks - 1)) {
                    bBlockWidth = bLastBlockWidth;
                }
                commBlockWidth = bz;
                for (int kb = 0; kb < commnBlocks; kb++) {
                    if (commLastBlockWidth > 0 && kb == (commnBlocks - 1)) {
                        commBlockWidth = commLastBlockWidth;
                    }

                    for (int i = ib * bz; i < (ib * bz) + aBlockHeight; i++) {
                        for (int j = jb * bz; j < (jb * bz) + bBlockWidth;
                             j++) {
                            for (int k = kb * bz;
                                 k < (kb * bz) + commBlockWidth; k++) {
                                if (A[i][k] != 0 && B[k][j] != 0) {
                                    C[i][j] += A[i][k] * B[k][j];
                                }
                            }
                        }
                    }
                }
            }
        }

        return C;
    }

    /**
     * Block matrix multiplication.
     *
     * @param A       - Matrix
     * @param B       - Matrix
     * @param aHeight - height of matrix A
     * @param bWidth  - width of matrix B
     * @param comm    - size of the common face
     * @return
     */
    public static double[][] matrixMultiply(
        double[][] A, double[][] B, int aHeight, int bWidth, int comm, int bz) {

        double[][] C = new double[aHeight][bWidth];

        int aHeightBlocks = aHeight / bz; // size = Height of A
        int aLastBlockHeight = aHeight - (aHeightBlocks * bz);
        if (aLastBlockHeight > 0) {
            aHeightBlocks++;
        }

        int bWidthBlocks = bWidth / bz; // size = Width of B
        int bLastBlockWidth = bWidth - (bWidthBlocks * bz);
        if (bLastBlockWidth > 0) {
            bWidthBlocks++;
        }

        int commnBlocks = comm / bz; // size = Width of A or Height of B
        int commLastBlockWidth = comm - (commnBlocks * bz);
        if (commLastBlockWidth > 0) {
            commnBlocks++;
        }

        int aBlockHeight = bz;
        int bBlockWidth = bz;
        int commBlockWidth = bz;

        for (int ib = 0; ib < aHeightBlocks; ib++) {
            if (aLastBlockHeight > 0 && ib == (aHeightBlocks - 1)) {
                aBlockHeight = aLastBlockHeight;
            }
            bBlockWidth = bz;
            commBlockWidth = bz;
            for (int jb = 0; jb < bWidthBlocks; jb++) {
                if (bLastBlockWidth > 0 && jb == (bWidthBlocks - 1)) {
                    bBlockWidth = bLastBlockWidth;
                }
                commBlockWidth = bz;
                for (int kb = 0; kb < commnBlocks; kb++) {
                    if (commLastBlockWidth > 0 && kb == (commnBlocks - 1)) {
                        commBlockWidth = commLastBlockWidth;
                    }

                    for (int i = ib * bz; i < (ib * bz) + aBlockHeight; i++) {
                        for (int j = jb * bz; j < (jb * bz) + bBlockWidth;
                             j++) {
                            for (int k = kb * bz;
                                 k < (kb * bz) + commBlockWidth; k++) {
                                if (A[i][k] != 0 && B[k][j] != 0) {
                                    C[i][j] += A[i][k] * B[k][j];
                                }
                            }
                        }
                    }
                }
            }
        }

        return C;
    }

    public static double[][] naiveMM(
        double[][] A, double[][] B, int aHeight, int bWidth, int comm) {

        double[][] C = new double[aHeight][bWidth];

        for (int i = 0; i < aHeight; i++) {

            for (int j = 0; j < bWidth; j++) {

                for (int k = 0; k < comm; k++) {
                    if (A[i][k] != 0 && B[k][j] != 0) {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
        }
        return C;
    }

    public static double[][] add(double[][] m1, double[][] m2) {

        int m1Len = m1.length;
        int m2Len = m2.length;

        if (m1Len == 0 || m2Len == 0 || m1Len != m2Len) {
            throw new RuntimeException("Incompatible Matrices");
        }
        else if (m1[0].length != m2[0].length || m1[0].length == 0 ||
                 m2[0].length == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }

        int vecLen = m1[0].length;
        double[][] res = new double[m1Len][vecLen];

        for (int i = 0; i < m1Len; i++) {
            for (int j = 0; j < vecLen; j++) {
                res[i][j] = m1[i][j] + m2[i][j];
            }
        }

        return res;
    }

    public static double[][] sub(double[][] m1, double[][] m2) {

        int m1Len = m1.length;
        int m2Len = m2.length;

        if (m1Len == 0 || m2Len == 0 || m1Len != m2Len) {
            throw new RuntimeException("Incompatible Matrices");
        }
        else if (m1[0].length != m2[0].length || m1[0].length == 0 ||
                 m2[0].length == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }

        int vecLen = m1[0].length;
        double[][] res = new double[m1Len][vecLen];

        for (int i = 0; i < m1Len; i++) {
            for (int j = 0; j < vecLen; j++) {
                res[i][j] = m1[i][j] - m2[i][j];
            }
        }

        return res;
    }

    // addition a double value to a matrix.
    public static double[][] add(double[][] m1, double s) {

        int m1Len = m1.length;

        if (m1Len == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }
        else if (m1[0].length == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }

        int vecLen = m1[0].length;
        double[][] res = new double[m1Len][vecLen];

        for (int i = 0; i < m1Len; i++) {
            for (int j = 0; j < vecLen; j++) {
                res[i][j] = m1[i][j] + s;
            }
        }
        return res;
    }

    // addition a double value to a matrix.
    public static double[][] sub(double[][] m1, double s) {

        int m1Len = m1.length;

        if (m1Len == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }
        else if (m1[0].length == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }

        int vecLen = m1[0].length;
        double[][] res = new double[m1Len][vecLen];

        for (int i = 0; i < m1Len; i++) {
            for (int j = 0; j < vecLen; j++) {
                res[i][j] = m1[i][j] - s;
            }
        }
        return res;
    }

    public static double[][] mul(double[][] m1, double[][] m2) {

        int m1Len = m1.length;
        int m2Len = m2.length;

        if (m1Len == 0 || m2Len == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }
        else if (m1[0].length == 0 || m2[0].length == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }
        else if (m1[0].length != m2.length) {
            throw new RuntimeException(
                "Incompatible Matrices for multiplication");
        }

        int m1Width = m1[0].length;
        int m2Width = m2[0].length;

        double[][] res = new double[m1Len][m2Width];

        double val = 0;
        for (int i = 0; i < m1Len; i++) {
            for (int j = 0; j < m2Width; j++) {
                val = 0;
                for (int k = 0; k < m1Width; k++) {
                    val += m1[i][k] * m2[k][j];
                }
                res[i][j] = val;
            }
        }
        return res;
    }

    // multiply a double value to a matrix.
    public static double[][] mult(double[][] m1, double s) {
        int m1Len = m1.length;

        if (m1Len == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }
        else if (m1[0].length == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }
        int vecLen = m1[0].length;
        double[][] res = new double[m1Len][vecLen];

        for (int i = 0; i < m1Len; i++) {
            for (int j = 0; j < vecLen; j++) {
                res[i][j] = m1[i][j] * s;
            }
        }

        return res;
    }

    public static double[][] copy(double[][] m1) {

        int m1Len = m1.length;

        if (m1Len == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }
        else if (m1[0].length == 0) {
            throw new RuntimeException("Incompatible Matrices");
        }
        int vecLen = m1[0].length;
        double[][] res = new double[m1Len][vecLen];
        for (int i = 0; i < m1Len; i++) {
            System.arraycopy(m1[i], 0, res[i], 0, vecLen);
        }
        return res;
    }

    public static void main(String[] args) {
        short[][] A = //{{0, 1, 1, 1},
            {{1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}};
        double[] Adiag = {4, 4, 4};
        double[][] B = {{1, 1}, {1, 1}, {1, 1}, {1, 1}};
        double[][] C = matrixMultiply(A, Adiag, B, 3, 2, 4, 2, 1);
        for (int i = 0; i < C.length; ++i) {
            for (int j = 0; j < C[i].length; ++j) {
                System.out.print(C[i][j] + " ");
            }
            System.out.println();
        }
    }
}

