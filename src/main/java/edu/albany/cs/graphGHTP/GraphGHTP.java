package edu.albany.cs.graphGHTP;

import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.Utils;
import edu.albany.cs.headApprox.PCSFHead;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.tailApprox.PCSFTail;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import java.util.*;

/**
 * This is GraphGHTP algorithm corresponding to Algorithm2 in our paper.
 * 
 * Please read our paper to find more details.
 *
 * @author baojian bzhou6@albany.edu
 */
public class GraphGHTP {

	/** graphSize, number of nodes in the input graph */
	private final int graphSize;
	/** edges in our graph, should notice that the graph should be connected. */
	private final ArrayList<Integer[]> edges;
	/** we use identity edge cost 1.0D in our algorithm. */
	private final ArrayList<Double> edgeCosts;
	/** the total sparsity s */
	private final int s;
	/** the maximum number of connected components formed by the forest F. */
	private final int g;
	/** budget B. */
	private final double B;
	/** we randomly select a single node to initialize x0. */
	private final boolean singleNodeInitial;
	/** counts info */
	private final double[] c;
	/** function that will be used in our algorithm. */
	private final Function function;
	/** parameter eta in our algorithm. */
	private double eta = 1.0D;
	/** only use for testing algorithm */
	private int verboseLevel = 0;
	/** connected components of current graph */
	private ConnectedComponents cc;
	/** final vector x that we get in our algorithm. */
	private final ArrayRealVector x;
	private final int[] trueSubGraph;
	private final boolean isNonTransportation;
	private final double epsilon = 1e-6;

	/** save function values in each iteration. */
	public ArrayList<Double> fValues;
	public int[] supportX;
	public HashSet<Integer> resultNodesTail = null;
	public double funcValueTail = 0.0D;
	public double runTime = 0.0D;
	/** number of iterations */
	public int iterNum;

	/**
	 * GraphModelIHTP algorithm
	 *
	 * @param edges
	 *            the edges of graph.
	 * @param edgeCosts
	 *            the cost of corresponding edges.
	 * @param c
	 *            costs or counts in data.
	 * @param s
	 *            sparsity
	 * @param g
	 *            number of cc
	 * @param B
	 *            the cost budget
	 * @param t
	 *            number of iterations
	 * @param singleNodeInitial
	 *            true : initialize x0 with single node.
	 * @param trueSubGraph
	 *            the true subgraph, in some data, there does not exist true
	 *            subgraph.
	 * @param func
	 *            used func.
	 * @param resultFileName
	 *            resultFileName in order to save file.
	 * @param fileName
	 *            fileName
	 */
	public GraphGHTP(int graphSize, ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, double[] c, int s, int g,
			double B, boolean singleNodeInitial, int[] trueSubGraph, Function func, boolean isNonTransPortation,
			double eta) {
		this.eta = eta;
		this.edges = edges;
		this.graphSize = graphSize;
		this.edgeCosts = edgeCosts;
		this.c = c;
		this.s = s;
		this.g = g;
		this.B = B;
		this.singleNodeInitial = singleNodeInitial;
		this.trueSubGraph = trueSubGraph;
		this.function = func;
		this.isNonTransportation = isNonTransPortation;
		if (checkInput()) {
			x = run(); // run the algorithm
		} else {
			x = null;
			System.out.println("input parameter is invalid.");
			System.exit(0);
		}
		if (verboseLevel >= 1) {
			printInfo();
		}
	}

	/**
	 * @return true if the input parameters are valid.
	 */
	private boolean checkInput() {
		Set<Integer> nodes = new HashSet<>();
		for (Integer[] edge : edges) {
			nodes.add(edge[0]);
			nodes.add(edge[1]);
		}
		/** nodes is inconsistent */
		if (nodes.size() != graphSize) {
			return false;
		}
		ArrayList<ArrayList<Integer>> adj = new ArrayList<>();
		for (int i = 0; i < graphSize; i++) {
			adj.add(new ArrayList<Integer>());
		}
		for (Integer[] edge : this.edges) {
			adj.get(edge[0]).add(edge[1]);
			adj.get(edge[1]).add(edge[0]);
		}
		cc = new ConnectedComponents(adj);
		return cc.connectivity;
	}

	private ArrayRealVector run() {
		long startTime = System.nanoTime();
		ArrayRealVector xi;
		if (singleNodeInitial) {
			/** return X0 with only one entry has 1.0D value. */
			xi = initializeXiRandomSingleNode();
		} else {
			/**
			 * return X0 with maximum connected component nodes have 1.0D
			 * values.
			 */
			xi = initializeXiMaximumCC(isNonTransportation);
			// xi = initializeXiMaximumCC();
		}

		fValues = new ArrayList<>();
		fValues.add(function.getFuncValue(xi.toArray()));
		while (true) {
			iterNum++;
			/** Gradient for the function. */
			double[] gradientFxi = function.getGradient(xi.toArray());
			gradientFxi = normalizeGradient(xi.toArray(), gradientFxi);
			/** head approximation */
			long startTimeHead = System.nanoTime();
			PCSFHead pcsfHead = new PCSFHead(edges, edgeCosts, gradientFxi, s, g, B, trueSubGraph);
			if (verboseLevel > 0) {
				System.out.println("Head runningTime: " + (System.nanoTime() - startTimeHead) / 1e9);
				System.out.println();
				for (int i : pcsfHead.bestForest.nodesInF) {
					System.out.print(gradientFxi[i] + " ");
				}
				System.out.println();
				for (int i = 0; i < 20; i++) {
					System.out.print(gradientFxi[i] + " ");
				}
				System.out.println();
				for (int i = 0; i < 20; i++) {
					System.out.print(gradientFxi[i] + " ");
				}
				System.out.println();
				PreRec preRec = new PreRec(pcsfHead.bestForest.nodesInF, trueSubGraph);
				System.out.println("size Head: " + pcsfHead.bestForest.nodesInF.size());
				System.out.println("trueSubGraph: " + trueSubGraph.length);
				System.out.println("pre: " + preRec.pre);
				System.out.println("rec: " + preRec.rec);
				// Utils.stop();
			}
			/** get head projection vector. */
			ArrayRealVector projectedHeader = projectVector(gradientFxi, pcsfHead.bestForest.nodesInF);
			/** get yi */
			ArrayRealVector yi = (new ArrayRealVector(xi)).add(projectedHeader.mapMultiply(eta));
			ArrayList<Integer> S = supp(yi);
			/** get arg min f(x) such that x_{S^c} = 0 */
			double[] b = function.getArgMaxFx(S);
			/** tail approximation */
			long startTimeTail = System.nanoTime();
			PCSFTail pcsfTail = new PCSFTail(edges, edgeCosts, b, s, g, B, trueSubGraph);
			xi = projectVector(b, pcsfTail.bestForest.nodesInF);
			resultNodesTail = new HashSet<>(pcsfTail.bestForest.nodesInF);
			double updatedFuncValue = getFuncTailVal(resultNodesTail);
			if (verboseLevel > 0) {
				System.out.println("updatedValue: " + updatedFuncValue);
				System.out.println("------------------------------------------------");
				System.out.println("Tail runningTime: " + (System.nanoTime() - startTimeTail) / 1e9);
				PreRec preRec = new PreRec(pcsfHead.bestForest.nodesInF, trueSubGraph);
				System.out.println("pre: " + preRec.pre);
				System.out.println("rec: " + preRec.rec);
				System.out.println("tail size: " + pcsfTail.bestForest.nodesInF.size());
				System.out.println("updatedFuncValue: " + updatedFuncValue);
				System.out.println("difference: " + (updatedFuncValue - fValues.get(fValues.size() - 1)));
				Utils.stop();
			}
			if ((updatedFuncValue - fValues.get(fValues.size() - 1)) < epsilon) {
				fValues.add(updatedFuncValue);
				break;
			} else {
				fValues.add(updatedFuncValue);
			}
		}
		supportX = getSupportNodes(xi.toArray());
		runTime = (System.nanoTime() - startTime) / 1e9;
		funcValueTail = getFuncTailVal(resultNodesTail);
		return xi;
	}

	private double[] normalizeGradient(double[] x, double[] gradient) {
		double[] normalizedGradient = new double[graphSize];
		for (int i = 0; i < graphSize; i++) {
			if ((gradient[i] < 0.0D) && (x[i] == 0.0D)) {
				normalizedGradient[i] = 0.0D;
			} else if ((gradient[i] > 0.0D) && (x[i] == 1.0D)) {
				normalizedGradient[i] = 0.0D;
			} else {
				normalizedGradient[i] = gradient[i];
			}
		}
		return normalizedGradient;
	}

	private double getFuncTailVal(HashSet<Integer> nodesInF) {
		double[] tmpX = new double[graphSize];
		Arrays.fill(tmpX, 0.0D);
		for (int i : nodesInF) {
			tmpX[i] = 1.0D;
		}
		return function.getFuncValue(tmpX);
	}

	private ArrayRealVector initializeXiMaximumCC(boolean isNonTransPortation) {
		int[] abnormalNodes = null;
		double mean = StatUtils.mean(c);
		double std = Math.sqrt(StatUtils.variance(c));
		if (isNonTransPortation) {
			for (int i = 0; i < graphSize; i++) {
				/** TODO to make sure mean + 2* std is for whole vector */
				if (Math.abs(c[i]) >= mean + 2.0D * std) {
					abnormalNodes = ArrayUtils.add(abnormalNodes, i);
				}
			}
		} else {
			for (int i = 0; i < graphSize; i++) {
				/** TODO to make sure mean + 2* std is for whole vector */
				if (Math.abs(c[i]) <= mean - 2.0D * std) {
					abnormalNodes = ArrayUtils.add(abnormalNodes, i);
				}
			}
		}
		if (abnormalNodes == null) {
			System.out.println("warning: the initial abnormal nodes is null ...");
			int maxIndex = 0;
			double maximalVal = -Double.MAX_VALUE;
			for (int i = 0; i < c.length; i++) {
				if (c[i] > maximalVal) {
					maximalVal = c[i];
					maxIndex = i;
				}
			}
			abnormalNodes = new int[] { maxIndex };
		}
		cc.computeCCSubGraph(abnormalNodes);
		int[] largestCC = cc.findLargestConnectedComponet(abnormalNodes);
		double[] x0 = new double[this.c.length];
		Arrays.fill(x0, 0.0D);
		for (int i = 0; i < largestCC.length; i++) {
			x0[largestCC[i]] = 1.0D;
		}
		if (verboseLevel > 0) {
			System.out.println("size of largestCC: " + largestCC.length);
			PreRec preRec = new PreRec(largestCC, trueSubGraph);
			System.out.println(preRec.toString());
		}
		return new ArrayRealVector(x0);
	}

	private ArrayRealVector initializeXiRandomSingleNode() {
		int[] abnormalNodes = null;
		double mean = StatUtils.mean(c);
		double std = Math.sqrt(StatUtils.variance(c));
		for (int i = 0; i < graphSize; i++) {
			if (c[i] >= mean + std) { // TODO need to change further
				abnormalNodes = ArrayUtils.add(abnormalNodes, i);
			}
		}

		double[] x0 = new double[c.length];
		Arrays.fill(x0, 0.0D);
		int index;
		if (abnormalNodes == null) {
			index = new Random().nextInt(x0.length);
			x0[index] = 1.0D;
		} else {
			index = new Random().nextInt(abnormalNodes.length);
			x0[abnormalNodes[index]] = 1.0D;
		}
		return new ArrayRealVector(x0);
	}

	/**
	 * get the support of a vector x
	 *
	 * @return a subset of nodes corresponding the index of vector x with
	 *         entries not equal to zero
	 */
	public ArrayList<Integer> supp(ArrayRealVector x) {
		if (x == null) {
			return null;
		}
		ArrayList<Integer> nodes = new ArrayList<Integer>();
		for (int i = 0; i < x.getDimension(); i++) {
			if (x.getEntry(i) != 0.0D) {
				nodes.add(i);
			}
		}
		return nodes;
	}

	/**
	 * get the nodes returned by algorithm
	 *
	 * @return the result nodes
	 */
	private int[] getSupportNodes(double[] x) {
		int[] nodes = null;
		for (int i = 0; i < x.length; i++) {
			if (x[i] != 0.0D) {
				nodes = ArrayUtils.add(nodes, i); // get nonzero nodes
			}
		}
		return nodes;
	}

	private ArrayRealVector projectVector(double[] x, ArrayList<Integer> projectionSet) {
		ArrayRealVector result = new ArrayRealVector(x);
		if (projectionSet == null) {
			Arrays.fill(x, 0.0D);
			return new ArrayRealVector(x);
		} else {
			for (int i = 0; i < result.getDimension(); i++) {
				if (!projectionSet.contains(i)) {
					result.setEntry(i, 0.0D);
				}
			}
		}
		return result;
	}

	public int[] getTailNodes() {
		int[] nodes = null;
		for (int k : this.resultNodesTail) {
			nodes = ArrayUtils.add(nodes, k);
		}
		return nodes;
	}

	public void printInfo() {
		System.out.println("total iterations : " + iterNum + "\n" + "final vector x : " + x.toString());
	}

	public static void main(String args[]) {
		ArrayRealVector arr = new ArrayRealVector(new double[] { 0.0, 1.0, 0.0 });
		for (int i = 0; i < arr.getDimension(); i++) {
			if (arr.getEntry(i) == 0.0D) {
				System.out.println("i : " + i + " " + arr.getEntry(i));
			}
			System.out.println(arr.toString());
		}
	}

}
