package edu.albany.cs.base;

import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.Edge;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * This class bean contains the input data
 * 
 * @author baojian bzhou6@albany.edu
 *
 */
public class InputData {

    public int numNodes;
    public int numEdges;
    public HashMap<Integer, Double> nodes;
    public HashMap<int[], Double> edges;
    public ArrayList<Integer[]> intEdges;
    public ArrayList<Double> edgeCosts;
    public ArrayList<Double> identityEdgeCosts;
    public ArrayList<Edge> newEdges;

    public String dataSource;
    public String usedAlgorithm;
    public int[] trueSubGraphNodes;
    public HashMap<int[], Double> trueSubGraphEdges = null;
    public ArrayList<ArrayList<Integer>> graphAdjList;
    public int[][] graphAdj;
    public ArrayList<ArrayList<Double>> graphWeightedAdjList;
    public double[][] graphWeightedAdj;
    public boolean connective;
    public int[] V;


    // for civilUnrest water network
    public double[] PValue;
    public double[] counts;
    // for citHepPh
    public double[] averCounts;
    // for transportation data
    public double[] speed;
    public double[] meanSpeed;
    public double[] speedStd;
    // for disease outbreak data
    public double[] mean;
    public double[] std;
    // connected components
    public ConnectedComponents cc;

    // identity base
    public double[] identityb;
    public double[] smallb;
}