package edu.albany.cs.headApprox;

/******************************************************************************
 * Compilation:  javac Digraph.java
 * Execution:    java Digraph filename.txt
 * Dependencies: Bag.java In.java StdOut.java
 * Data files:   http://algs4.cs.princeton.edu/42directed/tinyDG.txt
 * <p>
 * A graph, implemented using an array of lists.
 * Parallel edges and self-loops are permitted.
 * <p>
 * % java Digraph tinyDG.txt
 * 13 vertices, 22 edges
 * 0: 5 1
 * 1:
 * 2: 0 3
 * 3: 5 2
 * 4: 3 2
 * 5: 4
 * 6: 9 4 8 0
 * 7: 6 9
 * 8: 6
 * 9: 11 10
 * 10: 12
 * 11: 4 12
 * 12: 9
 ******************************************************************************/

import java.util.Stack;

/**
 *  The <tt>Digraph</tt> class represents a directed graph of vertices
 *  named 0 through <em>V</em> - 1.
 *  It supports the following two primary operations: add an edge to the digraph,
 *  iterate over all of the vertices adjacent from a given vertex.
 *  Parallel edges and self-loops are permitted.
 *  <p>
 *  This implementation uses an adjacency-lists representation, which 
 *  is a vertex-indexed array of {@link Bag} objects.
 *  All operations take constant time (in the worst case) except
 *  iterating over the vertices adjacent from a given vertex, which takes
 *  time proportional to the number of such vertices.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/42directed">Section 4.2</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */

public class Digraph {
    private static final String NEWLINE = System.getProperty("line.separator");

    private final int V;
    private int E;
    private Bag<Integer>[] adj;

    /**
     * Initializes an empty digraph with <em>V</em> vertices.
     *
     * @param  V the number of vertices
     * @throws IllegalArgumentException if V < 0
     */
    @SuppressWarnings("unchecked")
    public Digraph(int V) {
        if (V < 0) throw new IllegalArgumentException("Number of vertices in a Digraph must be nonnegative");
        this.V = V;
        this.E = 0;
        adj = (Bag<Integer>[]) new Bag[V];
        for (int v = 0; v < V; v++) {
            adj[v] = new Bag<Integer>();
        }
    }


    /**
     * Initializes a new digraph that is a deep copy of <tt>G</tt>.
     *
     * @param  G the digraph to copy
     */
    public Digraph(Digraph G) {
        this(G.V());
        this.E = G.E();
        for (int v = 0; v < G.V(); v++) {
            // reverse so that adjacency list is in same order as original
            Stack<Integer> reverse = new Stack<Integer>();
            for (int w : G.adj[v]) {
                reverse.push(w);
            }
            for (int w : reverse) {
                adj[v].add(w);
            }
        }
    }

    /**
     * Returns the number of vertices in this digraph.
     *
     * @return the number of vertices in this digraph
     */
    public int V() {
        return V;
    }

    /**
     * Returns the number of edges in this digraph.
     *
     * @return the number of edges in this digraph
     */
    public int E() {
        return E;
    }


    // throw an IndexOutOfBoundsException unless 0 <= v < V
    private void validateVertex(int v) {
        if (v < 0 || v >= V)
            throw new IndexOutOfBoundsException("vertex " + v + " is not between 0 and " + (V - 1));
    }

    /**
     * Adds the directed edge v->w to this digraph.
     *
     * @param  v the tail vertex
     * @param  w the head vertex
     * @throws IndexOutOfBoundsException unless both 0 <= v < V and 0 <= w < V
     */
    public void addEdge(int v, int w) {
        validateVertex(v);
        validateVertex(w);
        adj[v].add(w);
        E++;
    }

    /**
     * Returns the vertices adjacent from vertex <tt>v</tt> in this digraph.
     *
     * @param  v the vertex
     * @return the vertices adjacent from vertex <tt>v</tt> in this digraph, as an iterable
     * @throws IndexOutOfBoundsException unless 0 <= v < V
     */
    public Iterable<Integer> adj(int v) {
        validateVertex(v);
        return adj[v];
    }

    /**
     * Returns the number of directed edges incident from vertex <tt>v</tt>.
     * This is known as the <em>outdegree</em> of vertex <tt>v</tt>.
     *
     * @param  v the vertex
     * @return the outdegree of vertex <tt>v</tt>               
     * @throws IndexOutOfBoundsException unless 0 <= v < V
     */
    public int outdegree(int v) {
        validateVertex(v);
        return adj[v].size();
    }

    /**
     * Returns the reverse of the digraph.
     *
     * @return the reverse of the digraph
     */
    public Digraph reverse() {
        Digraph R = new Digraph(V);
        for (int v = 0; v < V; v++) {
            for (int w : adj(v)) {
                R.addEdge(w, v);
            }
        }
        return R;
    }

    /**
     * Returns a string representation of the graph.
     *
     * @return the number of vertices <em>V</em>, followed by the number of edges <em>E</em>,  
     *         followed by the <em>V</em> adjacency lists
     */
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(V + " vertices, " + E + " edges " + NEWLINE);
        for (int v = 0; v < V; v++) {
            s.append(String.format("%d: ", v));
            for (int w : adj[v]) {
                s.append(String.format("%d ", w));
            }
            s.append(NEWLINE);
        }
        return s.toString();
    }

    /**
     * Unit tests the <tt>Digraph</tt> data type.
     */
    public static void main(String[] args) {

    }

}
