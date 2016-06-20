package edu.albany.cs.fastPCST;

/**
 * Data pair
 *
 * @param <ValueType>
 * @param <IndexType>
 * @author baojian bzhou6@albany.edu
 */
public class Pair<ValueType, IndexType> {

    public ValueType value;
    public IndexType index;

    public Pair(ValueType value, IndexType index) {
        super();
        this.value = value;
        this.index = index;
    }

    public ValueType getFirst() {
        return value;
    }

    public IndexType getSecond() {
        return index;
    }

    @Override
    public String toString() {
        return "Pair [value=" + value + ", index=" + index + "]";
    }


}
