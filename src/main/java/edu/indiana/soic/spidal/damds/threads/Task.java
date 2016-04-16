package edu.indiana.soic.spidal.damds.threads;

public interface Task<T> {
    void run(T t);
}
