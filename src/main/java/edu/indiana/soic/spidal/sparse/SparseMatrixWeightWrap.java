package edu.indiana.soic.spidal.sparse;

public class SparseMatrixWeightWrap{

    SparseMatrix weight;
    SparseMatrix distance;
    private double avgDist = 1.0;
    private boolean isSammon;

    public SparseMatrixWeightWrap(SparseMatrix weight, SparseMatrix distance, boolean isSammon) {
        this.weight = weight;
        this.distance = distance;
        this.isSammon = isSammon;
    }

    public double getWeight(int columArrayIndex){
        double w = weight.getValues()[columArrayIndex];

        if(!isSammon) return w;

        double d = distance.getValues()[columArrayIndex];
        return w / Math.max(d, 0.001 * avgDist);
    }

    public void setAvgDistForSammon(double avgDist) {
        this.avgDist = avgDist;
    }

}
