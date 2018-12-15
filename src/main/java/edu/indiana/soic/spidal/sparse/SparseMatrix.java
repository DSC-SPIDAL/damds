package edu.indiana.soic.spidal.sparse;

/**
 * The sparse matrix is stored in CSR format
 * [Values], [Column indices] and [Row pointers]
 */
public class SparseMatrix {
    private final double[] values;
    private final int[] columns;
    private final int[] rowPointers;


    public SparseMatrix(double[] values, int[] columns, int[] rowPointers) {
        this.rowPointers = rowPointers;
        this.values = values;
        this.columns = columns;
    }

    public double[] getValues() {
        return values;
    }

    public int[] getColumns() {
        return columns;
    }

    public int[] getRowPointers() {
        return rowPointers;
    }
}