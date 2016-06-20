package edu.albany.cs.base;

import java.util.Arrays;

/**
 * Edge class, this class stands for an undirected graph edge
 *
 * @author Baojian bzhou6@albany.edu
 */
public class Edge {

    public final Integer i;                // end point i
    public final Integer j;                // end point j
    public final Integer ID;            // each edge has a ID
    public final Double cost;            // each edge has a cost
    public final Integer[] pair;


    /**
     * The default constructor will construct one edge which end points are i and j, cost is 0, ID is 0
     *
     * @param i
     * @param j
     */
    public Edge(int i, int j, Integer id, double cost) {
        this.i = i;
        this.j = j;
        this.ID = id;
        this.cost = cost;
        this.pair = new Integer[]{i, j};
    }

    /**
     * The default constructor will construct one edge which end points are i and j with cost and ID
     *
     * @param i    : end point i
     * @param j    : end point j
     * @param ID   : edge ID , we do not need to use it now.
     * @param cost :edge's cost
     */
    public Edge(int i, int j, int ID, double cost, Integer[] pair) {
        this.i = i;
        this.j = j;
        this.ID = ID;
        this.cost = cost;
        this.pair = pair;
    }

    @Override
    public String toString() {
        return "Edge [i=" + i + ", j=" + j + ", ID=" + ID + ", cost=" + cost + ", pair=" + Arrays.toString(pair) + "]";
    }

}
