package edu.albany.cs.base;

import org.apache.commons.lang3.ArrayUtils;

import java.util.*;

/**
 * ConnectedComponents is for computing the number of connected components of
 * graph.
 *
 * @author Baojian
 */
public class ConnectedComponents {

    // these private variables are used for computing the connected components of
    // subgraph induced by Integer[] A
    private Boolean[] marked = null; // marked[v] = has vertex v been marked?
    private Integer[] id = null; // id[v] = id of connected component containing v
    private Integer[] size = null; // size[id] = number of vertices in given
    // component
    private int count = 0; // number of connected components
    private ArrayList<ArrayList<Integer>> members = null; // this members is
    // connected components
    // nodes in subgraph
    private ArrayList<ArrayList<Integer>> graphAdjList;
    private HashMap<Integer, Integer> hashMap = null;

    public ArrayList<ArrayList<Integer>> components = null; // connected
    // components for
    // current subset
    public int numConnectedComponents = 0; // number of the connected components
    // in whole graph

    public boolean connectivity;

    public ArrayList<ArrayList<Integer>> updatedGraphAdj = new ArrayList<ArrayList<Integer>>();

    /**
     * given the adjacency list , we can calculate the connected components
     *
     * @param adjList the graph adjacency list
     */
    public ConnectedComponents(ArrayList<ArrayList<Integer>> adjList) {

        graphAdjList = adjList;
        numConnectedComponents = computeCCWholeGraph(graphAdjList);

        if (this.numConnectedComponents == 1) {
            connectivity = true;
        } else {
            connectivity = false;
        }
    }

    public boolean isConnected() {
        return connectivity;
    }

    /**
     * given the adjacency list , we can calculate the connected components
     *
     * @param adjList the graph adjacency list
     */
    public ConnectedComponents(int[][] adjList) {
        ArrayList<ArrayList<Integer>> arr = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < adjList.length; i++) {
            ArrayList<Integer> tt = new ArrayList<Integer>();
            for (int j : adjList[i]) {
                tt.add(j);
            }
            arr.add(tt);
        }
        this.graphAdjList = arr;
        this.numConnectedComponents = computeCCWholeGraph(this.graphAdjList);
        if (this.numConnectedComponents == 1) {
            connectivity = true;
        } else {
            connectivity = false;
        }
    }


    /**
     * given the adjacency matrix, we can calculate the connected components
     *
     * @param graph adjacency matrix
     */
    public ConnectedComponents(double[][] graph) {
        this.graphAdjList = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < graph.length; i++) {
            ArrayList<Integer> arr = new ArrayList<Integer>();
            for (int j = 0; j < graph.length; j++) {
                if (graph[i][j] >= 0) {
                    arr.add(j);
                }
            }
            this.graphAdjList.add(arr);
        }
        this.numConnectedComponents = computeCCWholeGraph(this.graphAdjList);

        if (this.numConnectedComponents == 1) {
            connectivity = true;
        } else {
            connectivity = false;
        }
    }

    /**
     * @param wholeGraphList the graph adjacency list
     * @return
     */
    public int computeCCWholeGraph(ArrayList<ArrayList<Integer>> wholeGraphList) {

        /**
         * these private variables are used for computing the connected components
         * of subgraph induced by Integer[] A
         */
        int numNodes = wholeGraphList.size();
        marked = new Boolean[numNodes]; // marked[v] = has vertex v been marked?
        count = 0; // number of connected components
        id = new Integer[numNodes];
        size = new Integer[numNodes];
        Arrays.fill(marked, false);
        Arrays.fill(id, 0);
        Arrays.fill(size, 0);

        int numOfCC = 0;
        for (int v = 0; v < numNodes; v++) {
            if (!marked[v]) {
                DFS(wholeGraphList, v);
                numOfCC++;
            }
        }
        // for the components
        Integer[] V = new Integer[numNodes];
        for (int i = 0; i < V.length; i++) {
            V[i] = i;
        }
        this.computeCCSubGraph(V);
        return numOfCC;
    }

    /**
     * check sub graph connectivity
     *
     * @param A the induced sub graph
     * @return if the induced sub graph is connected return true, else return
     * false
     */
    public boolean checkConnectivity(Integer[] A) {
        ArrayList<ArrayList<Integer>> subGraphList = generateSubGraph(A);
        Boolean[] marked = new Boolean[A.length]; // marked[v] = has vertex v been
        // marked?
        Integer[] id = new Integer[A.length]; // id[v] = id of connected component
        // containing v
        Integer[] size = new Integer[A.length]; // size[id] = number of vertices in
        // given component
        int count = 0; // number of connected components
        for (int v = 0; v < A.length; v++) {
            if (!marked[v]) {
                {
                    Stack<Integer> stack = new Stack<Integer>();
                    stack.push(v);
                    while (!stack.isEmpty()) {
                        v = (int) stack.pop();
                        if (marked[v] == false) {
                            marked[v] = true;
                            id[v] = count;
                            size[count]++;
                            for (int i : subGraphList.get(v)) {
                                stack.push(i);
                            }
                        }
                    }
                    count++; // we have added a new connected components
                    break;
                }
            }
        }

        if (A.length == size[0]) {
            return true;
        } else {
            return false;
        }
    }


    /**
     * check the this graph, if it is not connected then update the graph
     * structure by adding the random edges
     */
    public boolean checkConnectivity() {
        if (this.numConnectedComponents != 1) {
            System.out
                    .println("this graph is not connected, number of connected components in graph is : "
                            + this.numConnectedComponents);
            int num = this.graphAdjList.size();
            Integer[] V = new Integer[num];
            for (int i = 0; i < V.length; i++) {
                V[i] = i;
            }
            System.out.println("number of connected components is : " + this.computeCCSubGraph(V));
            // update the graph and let the graph to be connected
            int count = 0;
            ArrayList<Integer> largestConnecteComponent = this.components.get(0);
            for (ArrayList<Integer> arr : this.components) {
                if (count == 0) {
                    count++;
                    continue;
                }
                int root = largestConnecteComponent.get(count);
                int currentNeedToConnect = arr.get(0);
                this.graphAdjList.get(root).add(currentNeedToConnect);
                this.graphAdjList.get(currentNeedToConnect).add(root);
                count++;
            }
            System.out.println(
                    "after update the graph structure, the number of connected components in graph is : "
                            + this.computeCCSubGraph(V));
            return false;
        } else {
            return true;
        }
    }

    /**
     * check if this graph( edges ) is connected.
     *
     * @param edges represents a graph
     * @return
     */
    public final boolean checkConnectivity(ArrayList<Integer[]> edges) {
        Set<Integer> nodes = new HashSet<Integer>();
        for (Integer[] edge : edges) {
            nodes.add(edge[0]);
            nodes.add(edge[0]);
        }
        ArrayList<ArrayList<Integer>> adj = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < nodes.size(); i++) {
            adj.add(new ArrayList<Integer>());
        }
        for (Integer[] edge : edges) {
            adj.get(edge[0]).add(edge[1]);
            adj.get(edge[1]).add(edge[0]);
        }
        ConnectedComponents cc = new ConnectedComponents(adj);
        return cc.connectivity;
    }

    /**
     * @param A the reduced subset of nodes
     * @return the number of connected components in the reduced sub-graph by
     * subset A.
     */
    public int computeCCSubGraph(int[] A) {
        return this.computeCCSubGraph(Utils.intToInteger(A));
    }

    /**
     * @param A the reduced subset of nodes
     * @return the number of connected components in the reduced sub-graph by
     * subset A.
     */
    public int computeCCSubGraph(ArrayList<Integer> A) {
        Integer[] nodes = new Integer[A.size()];
        for (int i = 0; i < A.size(); i++) {
            nodes[i] = A.get(i);
        }
        return this.computeCCSubGraph(nodes);
    }

    /**
     * @param A the reduced subset of nodes
     * @return the number of connected components in the reduced sub-graph by
     * subset A.
     */
    public int computeCCSubGraph(Integer[] A) {
        if (A == null) {
            return 0;
        }
        // these private variables are used for computing the connected components
        // of subgraph induced by Integer[] A
        marked = null; // marked[v] = has vertex v been marked?
        id = null; // id[v] = id of connected component containing v
        size = null; // size[id] = number of vertices in given component
        count = 0; // number of connected components

        ArrayList<ArrayList<Integer>> subGraphList = generateSubGraph(A);
        hashMap = new HashMap<Integer, Integer>();
        for (int k = 0; k < A.length; k++) {
            hashMap.put(k, A[k]);
        }
        int numNodes = subGraphList.size();
        marked = new Boolean[numNodes];
        id = new Integer[numNodes];
        size = new Integer[numNodes];

        for (int k = 0; k < marked.length; k++) {
            marked[k] = false;
        }
        for (int k = 0; k < id.length; k++) {
            id[k] = 0;
        }
        for (int k = 0; k < size.length; k++) {
            size[k] = 0;
        }

        for (int v = 0; v < numNodes; v++) {
            if (!marked[v]) {
                DFS(subGraphList, v);
                count++; // we have added a new connected components
            }
        }

        {// this task is to find the members of each connected components
            members = new ArrayList<ArrayList<Integer>>(count);
            for (int i = 0; i < count; i++) { // compute members
                ArrayList<Integer> temp = new ArrayList<Integer>();
                for (int j = 0; j < id.length; j++) {
                    if (id[j] == i) {
                        temp.add(j);
                    }
                }
                members.add(temp);
            }

            ArrayList<ArrayList<Integer>> comps = new ArrayList<ArrayList<Integer>>(count); // map
            // to
            // our
            // original
            // nodes
            // ID
            for (int i = 0; i < count; i++) {
                ArrayList<Integer> temp = new ArrayList<Integer>();
                for (int k : members.get(i)) {
                    temp.add(hashMap.get(k));
                }
                comps.add(temp);
            }
            this.components = comps;
        }
        return count;
    }


    /**
     * @return the largest connected components in current graph
     */
    public int[] findLargestConnectedComponet(int[] A) {
        this.computeCCSubGraph(A);
        if (this.components != null) {
            int count = Integer.MIN_VALUE;
            int[] largest = null;
            for (ArrayList<Integer> mem : this.components) {
                if (count < mem.size()) {
                    count = mem.size();
                    largest = null;
                    for (int i : mem) {
                        largest = ArrayUtils.add(largest, i);
                    }
                }
            }
            return largest;
        } else {
            return null;
        }
    }

    /**
     * @return the largest connected components in current graph
     */
    public int[] findLargestConnectedComponet(ArrayList<Integer> A) {
        int[] AA = new int[A.size()];
        for (int i = 0; i < A.size(); i++) {
            AA[i] = A.get(i);
        }
        this.computeCCSubGraph(AA);
        if (this.components != null) {
            int count = Integer.MIN_VALUE;
            int[] largest = null;
            for (ArrayList<Integer> mem : this.components) {
                if (count < mem.size()) {
                    count = mem.size();
                    largest = null;
                    for (int i : mem) {
                        largest = ArrayUtils.add(largest, i);
                    }
                }
            }
            return largest;
        } else {
            return null;
        }
    }

    /**
     * according to the sub-graph of nodes, we create a subgraph of the graph G
     */
    private ArrayList<ArrayList<Integer>> generateSubGraph(Integer[] A) {
        if (A == null) {
            return null;
        }
        HashMap<Integer, Integer> hashMap = new HashMap<Integer, Integer>();
        for (int i = 0; i < A.length; i++) {
            hashMap.put(A[i], i);
        }
        ArrayList<ArrayList<Integer>> subGraphList = new ArrayList<ArrayList<Integer>>();
        for (int i : A) {
            ArrayList<Integer> temp = new ArrayList<Integer>(this.graphAdjList.get(i));
            temp.retainAll(Arrays.asList(A));
            subGraphList.add(temp);
        }
        // replace with the new sequences of the nodes
        for (ArrayList<Integer> arrayList : subGraphList) {
            for (int i = 0; i < arrayList.size(); i++) {
                arrayList.set(i, hashMap.get(arrayList.get(i)));
            }
        }
        return subGraphList;
    }

    /**
     * @param subGraphList reduced sub-graph adjacency list
     * @param v            starting node
     */
    private void DFS(ArrayList<ArrayList<Integer>> subGraphList, int v) {
        Stack<Integer> stack = new Stack<Integer>();
        stack.push(v);
        while (!stack.isEmpty()) {
            v = (int) stack.pop();
            if (marked[v] == false) {
                marked[v] = true;
                id[v] = count;
                size[count]++;
                for (int i : subGraphList.get(v)) {
                    stack.push(i);
                }
            }
        }
    }

    public static void main(String args[]) {
        APDMInputFormat apdm =
                new APDMInputFormat("simuDataSet/GridDataSet/Grid-Data-100/APDM-GridData-100-0.txt");
        // APDMInputFormat apdm = new
        // APDMInputFormat("realDataSet/Disease/Chile_2_year/2013-10-18.txt") ;
        ConnectedComponents c = new ConnectedComponents(apdm.data.graphAdjList);
        System.out.println("the number of connected components is : " + c.numConnectedComponents);
        c.computeCCSubGraph(Utils.intToInteger(apdm.data.V));
        System.out.println("the number of connected components is : " + c.numConnectedComponents);
        int count = 0;
        for (ArrayList<Integer> arr : c.components) {
            System.out.println(count++ + " : " + arr.toString());
        }
        System.out.println("[0] should be [1]: " + c.computeCCSubGraph(new Integer[]{0}));
        System.out.println("[1] should be [1]: " + c.computeCCSubGraph(new Integer[]{1}));
        System.out.println("[0 1] should be [1]: " + c.computeCCSubGraph(new Integer[]{0, 1}));
        System.out
                .println("[0 1 4 6] should be [3]: " + c.computeCCSubGraph(new Integer[]{0, 1, 4, 6}));
        System.out.println(
                "[0 1 2 3 4 6] should be [2]: " + c.computeCCSubGraph(new Integer[]{0, 1, 2, 3, 4, 6}));
        System.out.println("[0 3 4 5 88 91 ] should be [4]: "
                + c.computeCCSubGraph(new Integer[]{0, 3, 4, 5, 88, 91}));
        for (ArrayList<Integer> com : c.components) {
            System.out.println(com.toString());
        }
    }
}
