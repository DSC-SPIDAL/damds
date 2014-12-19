package edu.indiana.soic.spidal.common;

import java.io.Serializable;

public final class Range implements Serializable {
    /**
     * The rangeStart index
     */
    private int startIndex;
    /**
     * The rangeEnd index.
     */
    private int endIndex;
    /**
     * The length of the range.
     */
    private int length;
    private String privateStartSeqName;
    private String privateEndSeqName;

    /**
     * Initializes a new instance of the class.
     *
     * @param start The starting index of the Range.
     * @param end   The ending index of the Range.
     */
    public Range(int start, int end) {
        startIndex = start;
        endIndex = end;
        length = end - start + 1;
    }

    public int getStartIndex() {
        return startIndex;
    }

    public int getEndIndex() {
        return endIndex;
    }

    public int getLength() {
        return length;
    }

    public String getStartSeqName() {
        return privateStartSeqName;
    }

    public void setStartSeqName(String value) {
        privateStartSeqName = value;
    }

    public String getEndSeqName() {
        return privateEndSeqName;
    }

    public void setEndSeqName(String value) {
        privateEndSeqName = value;
    }

    public boolean Contains(int index) {
        return (index >= startIndex && index <= endIndex);
    }

    /**
     * Returns the fully qualified type name of this instance.
     *
     * @return A <see cref="T:System.String"/> containing a fully qualified type name.
     */
    @Override
    public String toString() {
        return String.format("(%1$s:%2$s)", Integer.toString(startIndex), Integer.toString(endIndex));
    }

    /**
     * Returns true if there's an intersection of this range with the given range
     *
     * @param range The range to see if an intersection exists
     * @return True if an intersection exists or false otherwise
     */
    public boolean IntersectsWith(Range range) {
        Range lengthiest = range.getLength() >= length ? range : this;
        Range other = range == lengthiest ? this : range;

        return lengthiest.Contains(other.getStartIndex()) || lengthiest.Contains(other.getEndIndex());
    }

    /**
     * Gets the intersecting range assuming an intersection exists. Use <code>IntersectsWith(Range range)</code>
     * to check for an existing intersection
     *
     * @param range The range to intersect with
     * @return The intersection with the given range
     */
    public Range GetIntersectionWith(Range range) {
        int start = range.getStartIndex() >= startIndex ? range.getStartIndex() : getStartIndex();
        int end = range.getEndIndex() <= endIndex ? range.getEndIndex() : endIndex;
        if (start <= end) {
            return new Range(start, end);
        }
        throw new RuntimeException(String.format("Given range [%1$s, %2$s] does not intersect with this range [%3$s, " +
                "%4$s]", range.getStartIndex(), range.getEndIndex(), startIndex, endIndex));
    }
}

