package edu.indiana.soic.spidal.sparse;

public class SparseMatrixWeightWrap{

    SparseMatrix weight;
    SparseMatrix distance;
    private boolean isSammon;
    public SparseMatrixWeightWrap(SparseMatrix weight, SparseMatrix distance, boolean isSammon) {
        this.weight = weight;
        this.distance = distance;
        this.isSammon = isSammon;
    }

    public double getWeight(int i, int j){
        double w = 0;
        // get
        if(!isSammon) return w;

//        double d = distances[i][j] * 1.0 / Short.MAX_VALUE;
//        return w / Math.max(d, 0.001 * avgDist);
        return w;
    }
}
