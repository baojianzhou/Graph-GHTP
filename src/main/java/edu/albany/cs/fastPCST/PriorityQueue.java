package edu.albany.cs.fastPCST;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Priority Queue data structure
 * 
 * we set the value type to <Double,Integer>
 *
 * @author baojian bzhou6@albany.edu
 * 
 * 
 */
public class PriorityQueue {

	private SortedSet<Pair<Double, Integer>> sorted_set = new TreeSet<Pair<Double, Integer>>(new Comp());
	private ArrayList<Pair<Double, Integer>> index_to_iterator = new ArrayList<Pair<Double, Integer>>();

	// reference parameters
	public Double get_min_firstP = -1d;
	public Integer get_min_secondP = -1;
	public Double delete_min_firstP = -1d;
	public Integer delete_min_secondP = -1;

	public PriorityQueue() {
	}

	public boolean is_empty() {
		return sorted_set.isEmpty();
	}

	public boolean get_min(Double value, Integer index) {
		if (sorted_set == null || sorted_set.isEmpty()) {
			return false;
		}
		value = sorted_set.first().getFirst();
		index = sorted_set.first().getSecond();
		this.get_min_firstP = value;
		this.get_min_secondP = index;
		return true;
	}

	public boolean delete_min(Double value, Integer index) {
		if (sorted_set.isEmpty()) {
			return false;
		}
		value = sorted_set.first().getFirst();
		index = sorted_set.first().getSecond();
		sorted_set.remove(sorted_set.first());
		this.delete_min_firstP = value;
		this.delete_min_secondP = index;
		return true;
	}

	public void insert(Double value, Integer index) {
		if (index >= index_to_iterator.size()) {
			// resize the index_to_iterator
			int newSize = (index + 1) - index_to_iterator.size();
			if (index_to_iterator.size() > (index + 1)) {
				for (int i = 0; i < Math.abs(newSize); i++) {
					index_to_iterator.remove(index_to_iterator.size() - 1);
				}
			} else {
				for (int i = 0; i < newSize; i++) {
					index_to_iterator.add(new Pair<Double, Integer>(0.0, 1));
				}
			}
		}
		Pair<Double, Integer> pair = new Pair<Double, Integer>(value, index);
		sorted_set.add(pair);
		index_to_iterator.set(index, pair);
	}

	public final void decrease_key(Double new_value, Integer index) {
		sorted_set.remove(index_to_iterator.get(index));
		Pair<Double, Integer> pair = new Pair<Double, Integer>(new_value, index);
		sorted_set.add(pair);
		index_to_iterator.set(index, pair);
	}

	public final void delete_element(Integer index) {
		sorted_set.remove(index_to_iterator.get(index));
	}

	public int getQueueSize() {
		return sorted_set.size();
	}

	public void print() {
		int ii = 0;
		for (Pair<Double, Integer> pair : sorted_set) {
			System.out.format("(ii:%d) %f %d\n", ii, pair.getFirst(), pair.getSecond());
			ii++;
		}
	}

	public static class Comp implements Comparator<Pair<Double, Integer>> {
		// @Override
		public int compare(Pair<Double, Integer> pair1, Pair<Double, Integer> pair2) {
			if (pair1.equals(null) || pair2.equals(null)) {
				new IllegalArgumentException("cannot be null");
				System.exit(0);
				return -1;
			}
			if ((pair1.getFirst() - pair2.getFirst()) > 0) {
				return 1;
			} else if ((pair1.getFirst() - pair2.getFirst()) < 0) {
				return -1;
			} else {
				if ((pair1.getSecond() - pair2.getSecond()) > 0) {
					return 1;
				} else if ((pair1.getSecond() - pair2.getSecond()) < 0) {
					return -1;
				} else {
					return 0;
				}
			}
		}
	}

	// test
	public static void main(String args[]) {
		SortedSet<Pair<Double, Integer>> sorted_set = new TreeSet<Pair<Double, Integer>>(new Comp());
		sorted_set.add(new Pair<Double, Integer>(1.0, 1));
		sorted_set.add(new Pair<Double, Integer>(1.0, 2));
		sorted_set.add(new Pair<Double, Integer>(1.5, 1));
		int i = 0;
		for (Pair<Double, Integer> pair : sorted_set) {
			System.out.println("(i:" + i++ + ") " + pair.getFirst() + " " + pair.getSecond());
		}
	}
}