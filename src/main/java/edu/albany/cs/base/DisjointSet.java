package edu.albany.cs.base;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/*******************************************************************************
 * Copyright (c) 2015 IBM Corporation and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 * <p>
 * Contributors:
 * IBM Corporation - initial API and implementation
 * Andrey Loskutov <loskutov@gmx.de> - generified interface, bug 462760
 ******************************************************************************/


/**
 * A disjoint set is a generic data structure that represents a collection of
 * sets that are assumed to be disjoint (no object exists in more than one set).
 * <p>
 * This disjoint set implementation represents the disjoint set as a forest,
 * where the nodes of each tree all belong to the same set. This implementation
 * uses path compression in the findSet implementation to flatten each tree to a
 * constant depth. A rank is maintained for each tree that is used when
 * performing union operations to ensure the tree remains balanced.
 * <p>
 * Ref: Cormen, Leiserson, and Rivest <it>Introduction to Algorithms</it>,
 * McGraw-Hill, 1990. The disjoint set forest implementation in section 22.3.
 * </p>
 *
 * @param <T>
 *            element type
 * @since 3.2
 */
public class DisjointSet<T> {

	public int numConnectedComponents = 0;

	/**
	 * A node in the disjoint set forest. Each tree in the forest is a disjoint
	 * set, where the root of the tree is the set representative.
	 */
	private static class Node<T> {
		/** The node rank used for union by rank optimization */
		int rank;
		/** The parent of this node in the tree. */
		T parent;

		Node(T parent, int rank) {
			this.parent = parent;
			this.rank = rank;
		}
	}

	/**
	 * Map of Object -> Node, where each key is an object in the disjoint set,
	 * and the Node represents its position and rank within the set.
	 */
	private final HashMap<T, Node<T>> objectsToNodes = new HashMap<T, Node<T>>();

	/**
	 * Returns the set token for the given object, or null if the object does
	 * not belong to any set. All object in the same set have an identical set
	 * token.
	 * 
	 * @param o
	 *            The object to return the set token for
	 * @return The set token, or <code>null</code>
	 */
	public T findSet(Object o) {
		DisjointSet.Node<T> node = objectsToNodes.get(o);
		if (node == null) {
			return null;
		}
		if (o != node.parent) {
			node.parent = findSet(node.parent);
		}
		return node.parent;
	}

	/**
	 * added by Baojian
	 */
	public int numConnectedComponents() {
		return this.numConnectedComponents;
	}

	/**
	 * added by Baojian
	 */
	public HashMap<T, Set<T>> getConnectedComponents() {
		HashMap<T, Set<T>> hashMap = new HashMap<T, Set<T>>();
		for (T t : this.objectsToNodes.keySet()) {
			T t1 = findSet(t);
			if (!hashMap.containsKey(t1)) {
				Set<T> set = new HashSet<T>();
				set.add(t);
				hashMap.put(t1, set);
			} else {
				hashMap.get(t1).add(t);
			}
		}
		return hashMap;
	}

	/**
	 * added by Baojian
	 */
	public boolean isConnected(Object q, Object p) {
		if (q == null && p == null) {
			return false;
		}
		return findSet(q) == findSet(p);
	}

	/**
	 * Adds a new set to the group of disjoint sets for the given object. It is
	 * assumed that the object does not yet belong to any set.
	 * 
	 * @param o
	 *            The object to add to the set
	 */
	public void makeSet(T o) {
		objectsToNodes.put(o, new Node<T>(o, 0));
		numConnectedComponents++;
	}

	/**
	 * Removes all elements belonging to the set of the given object.
	 * 
	 * @param o
	 *            The object to remove
	 */
	public void removeSet(Object o) {
		Object set = findSet(o);
		if (set == null) {
			return;
		}
		for (Iterator<T> it = objectsToNodes.keySet().iterator(); it.hasNext();) {
			T next = it.next();
			// remove the set representative last, otherwise findSet will fail
			if (next != set && findSet(next) == set) {
				it.remove();
			}
		}
		objectsToNodes.remove(set);
	}

	/**
	 * Copies all objects in the disjoint set to the provided list
	 * 
	 * @param list
	 *            The list to copy objects into
	 */
	public void toList(List<? super T> list) {
		list.addAll(objectsToNodes.keySet());
	}

	/**
	 * Unions the set represented by token x with the set represented by token
	 * y. Has no effect if either x or y is not in the disjoint set, or if they
	 * already belong to the same set.
	 * 
	 * @param x
	 *            The first set to union
	 * @param y
	 *            The second set to union
	 */
	public void union(T x, T y) {
		T setX = findSet(x);
		T setY = findSet(y);
		if (setX == null || setY == null || setX == setY || setX.equals(setY)) {
			return;
		}
		Node<T> nodeX = objectsToNodes.get(setX);
		Node<T> nodeY = objectsToNodes.get(setY);
		// join the two sets by pointing the root of one at the root of the
		// other
		if (nodeX.rank > nodeY.rank) {
			nodeY.parent = x;
		} else {
			nodeX.parent = y;
			if (nodeX.rank == nodeY.rank) {
				nodeY.rank++;
			}
		}
		numConnectedComponents--;
	}

	public static void testCase1() {
		DisjointSet<Integer> dis = new DisjointSet<Integer>(); 
		dis.makeSet(0);
		dis.makeSet(1);
		dis.makeSet(2);
		dis.makeSet(3);
		dis.makeSet(4);
		dis.makeSet(5);
		dis.union(0, 1);
		dis.union(1, 2);
		dis.union(0, 2);
		dis.union(2, 0);
		dis.union(3, 4);
		for (Integer s : dis.getConnectedComponents().keySet()) {
			System.out.println("ID : " + s + " ; nodes : " + dis.getConnectedComponents().get(s));
		}
		System.out.println("number of connected components : " + dis.numConnectedComponents);
	}

	public static void testCase2() {
		DisjointSet<Integer> dis = new DisjointSet<Integer>(); 
		dis.makeSet(0);
		dis.makeSet(1);
		dis.makeSet(2);
		dis.makeSet(3);
		dis.makeSet(4);
		dis.makeSet(5);
		dis.union(0, 1);
		dis.union(1, 2);
		dis.union(0, 2);
		dis.union(2, 0);
		dis.union(3, 4);
		dis.union(3, 1);
		dis.union(3, 5);
		for (Integer s : dis.getConnectedComponents().keySet()) {
			System.out.println("ID : " + s + " ; nodes : " + dis.getConnectedComponents().get(s));
		}
		System.out.println("number of connected components : " + dis.numConnectedComponents);
	}

	public static void main(String args[]) {
		testCase1();
		testCase2();
	}

}