package edu.albany.cs.base;

import java.io.IOException;
import java.util.HashSet;

import org.apache.commons.lang3.ArrayUtils;

public final class Utils {

	private Utils() {

	}

	public final static String commentLine = "#################################################################";

	/**
	 * transform int array to Integer array
	 * 
	 * @param array
	 *            input array
	 * @return Integer array
	 */
	public static Integer[] intToInteger(int[] array) {
		if (array == null) {
			return null;
		} else {
			Integer[] newArray = new Integer[array.length];
			for (int i = 0; i < newArray.length; i++) {
				newArray[i] = array[i];
			}
			return newArray;
		}
	}

	/**
	 * 
	 * return the intersection of a and b
	 * @param a array of elements (may not be unify)
	 * @param b array of elements (may not be unify)
	 * @return a intersect b
	 */
	public static int[] intersect(int[] a, int[] b) {
		if (a == null || b == null) {
			return null;
		}
		int[] result = new int[Math.min(a.length, b.length)];
		HashSet<Integer> hs = new HashSet<Integer>();
		for (int i = 0; i < b.length; i++) {
			hs.add(b[i]);
		}
		int count = 0;
		for (int i = 0; i < a.length; i++) {
			if (hs.contains(a[i])) {
				result[count++] = a[i];
			}
		}
		return ArrayUtils.subarray(result, 0, count);
	}
	
	/**
	 * stop for testing
	 */
	public static void stop() {
		try {
			System.out.println("Press any key to continue ...");
			System.in.read();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static int[] setToArray(HashSet<Integer> set) {
		int[] result = null;
		for (int i : set) {
			result = ArrayUtils.add(result, i);
		}
		return result;
	}

}
