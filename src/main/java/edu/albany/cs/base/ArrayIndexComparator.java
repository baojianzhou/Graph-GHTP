package edu.albany.cs.base;

import org.apache.commons.lang3.ArrayUtils;

import java.util.Arrays;
import java.util.Comparator;

/**
 * @author baojian bzhou6@albany.edu
 *
 */
public class ArrayIndexComparator implements Comparator<Integer> {

	/** we want to sort array by descending order */
	private final Double[] array;
	/** keep indexes in descending order */
	public final Integer[] indexes;

	public ArrayIndexComparator(Double[] array) {
		this.array = array;
		this.indexes = createIndexArray();
	}

	public ArrayIndexComparator(double[] array) {
		this.array = ArrayUtils.toObject(array);
		this.indexes = createIndexArray();
	}

	public Integer[] createIndexArray() {
		Integer[] indexes = new Integer[array.length];
		for (int i = 0; i < this.array.length; i++) {
			indexes[i] = i;
		}
		return indexes;
	}

	/**
	 * keep the descending order i1>i2> ...
	 * 
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	@Override
	public int compare(Integer index1, Integer index2) {
		return array[index2].compareTo(array[index1]);
	}

	/**
	 * simple test
	 * @param args
	 */
	public static void main(String args[]) {
		Double[] vectorRatioCB = new Double[] { 0.1, 0.03, 0.2, 0.04, 0.5 };
		ArrayIndexComparator arrayIndexComparator = new ArrayIndexComparator(vectorRatioCB);
		Integer[] indexes = arrayIndexComparator.indexes;
		Arrays.sort(indexes, arrayIndexComparator);
		System.out.println(Arrays.toString(indexes));
		for (int i : indexes) {
			System.out.print(vectorRatioCB[i] + " ");
		}
	}
}