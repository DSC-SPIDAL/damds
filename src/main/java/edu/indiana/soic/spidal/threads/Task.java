package edu.indiana.soic.spidal.threads;

public interface Task<T> {
    void run(T t);
}
