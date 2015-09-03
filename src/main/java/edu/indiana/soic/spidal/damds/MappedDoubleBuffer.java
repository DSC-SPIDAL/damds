package edu.indiana.soic.spidal.damds;

import java.nio.DoubleBuffer;
import java.nio.MappedByteBuffer;

public class MappedDoubleBuffer{
    private MappedByteBuffer mbb;
    private DoubleBuffer db;

    public MappedDoubleBuffer(MappedByteBuffer mbb) {
        this.mbb = mbb;
        db = mbb.asDoubleBuffer();
    }

    public DoubleBuffer getDb() {
        return db;
    }

    public void force(){
        mbb.force();
    }

    public void position(int pos) {
        mbb.position(pos);
        db.position(pos);
    }
}
