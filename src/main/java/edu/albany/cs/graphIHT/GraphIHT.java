package edu.albany.cs.graphIHT;

import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.headApprox.PCSFHead;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.tailApprox.PCSFTail;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import java.util.*;

/**
 * Algorithm 1 : Graph-IHT algorithm
 *
 * @author Baojian Zhou (bzhou6@albany.edu)
 */
public class GraphIHT {

	/** graphSize, number of nodes in the input graph */
	private final int graphSize;
	/** edges in our graph, should notice that the graph should be connected. */
	private final ArrayList<Integer[]> edges;
	/** we use identity edge cost 1.0D in our algorithm. */
	private final ArrayList<Double> edgeCosts;
	/** total sparsity of S */
	private final int s;
	/** # of connected components in forest F */
	private final int g;
	/** bound on the total weight w(F) of edges in the */
	private final double B;
	/** we randomly select a single node to intialize x0 if it is true. */
	private final boolean singleNodeInitial;
	/** data */
	private final double[] c;
	/** statistic function */
	private final Function function;
	/** parameter delta in our algorithm */
	private final double eta = 1.0D;
	/** connected components of current graph */
	private ConnectedComponents cc;
	/** halting condition */
	private double epsilon = 1e-6;

	private final int[] trueSubGraph;
	private int verboseLevel = 0;

	public int[] supportX;
	/** final vector x that we get in our algorithm. */
	public ArrayRealVector x;
	public HashSet<Integer> resultNodesTail;
	public double funcValueTail = 0.0D;
	public double runTime = 0.0D;
	public ArrayList<Double> fValues;
	/** number of iterations */
	public int iterNum;

	public GraphIHT(int graphSize, ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, double[] c, int s, int g,
			double B, boolean singleNodeInitial, int[] trueSubGraph, Function func) {
		this.graphSize = graphSize;
		this.edgeCosts = edgeCosts;
		this.edges = edges;
		this.c = c;
		this.s = s;
		this.g = g;
		this.B = B;
		this.singleNodeInitial = singleNodeInitial;
		this.trueSubGraph = trueSubGraph;
		this.function = func;
		if (checkInput()) {
			x = run();
		} else {
			x = null;
			System.out.println("input parameter is invalid. ");
			System.exit(0);
		}
	}

	private ArrayRealVector run() {
		long startTime = System.nanoTime();
		ArrayRealVector xi;
		if (singleNodeInitial) {
			/** return X0 with only one entry has 1.0D value */
			xi = initializeX_RandomSingleNode();
		} else {
			/** return X0 with maximum connected component have 1.0D values. */
			xi = initializeX_MaximumCC();
		}

		fValues = new ArrayList<Double>();
		fValues.add(function.getFuncValue(xi.toArray()));
		while (true) {
			iterNum++;
			/** Gradient for the function. */
			double[] gradientFxi = function.getGradient(xi.toArray());
			/** TODO check this further */
			gradientFxi = normalizeGradient(xi.toArray(), gradientFxi);
			/** head approximation */
			PCSFHead pcsfHead = new PCSFHead(edges, edgeCosts, gradientFxi, s, g, B, trueSubGraph);
			/** get head projection vector. */
			ArrayRealVector projectedHeader = projectVector(gradientFxi, pcsfHead.bestForest.nodesInF);
			/** get yi */
			ArrayRealVector yi = (new ArrayRealVector(xi)).add(projectedHeader.mapMultiply(eta));
			/** tail approximation */
			PCSFTail pcsfTail = new PCSFTail(edges, edgeCosts, yi.toArray(), s, g, B, trueSubGraph);
			/** get x_{i+1} */
			xi = projectVector(yi.toArray(), pcsfTail.bestForest.nodesInF);
			xi = normalize(xi);
			resultNodesTail = new HashSet<>(pcsfTail.bestForest.nodesInF);
			double updatedFuncValue = getFuncTailVal(resultNodesTail);
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

	private ArrayRealVector initializeX_MaximumCC() {
		int[] abnormalNodes = null;
		double mean = StatUtils.mean(c);
		double std = Math.sqrt(StatUtils.variance(c));
		for (int i = 0; i < graphSize; i++) {
			if (c[i] >= mean + 2.0D * std) {
				abnormalNodes = ArrayUtils.add(abnormalNodes, i);
			}
		}
		cc.computeCCSubGraph(abnormalNodes);
		int[] largestCC = cc.findLargestConnectedComponet(abnormalNodes);
		double[] x0 = new double[graphSize];
		Arrays.fill(x0, 0.0D);
		for (int i = 0; i < largestCC.length; i++) {
			x0[largestCC[i]] = 1.0D;
		}
		if (verboseLevel > 0) {
			PreRec preRec = new PreRec(largestCC, trueSubGraph);
			System.out.println(preRec.toString());
		}
		return new ArrayRealVector(x0);
	}

	private ArrayRealVector initializeX_RandomSingleNode() {
		int[] abnormalNodes = null;
		double mean = StatUtils.mean(c);
		double std = Math.sqrt(StatUtils.variance(c));
		for (int i = 0; i < graphSize; i++) {
			if (c[i] >= mean + 2.0D * std) {
				abnormalNodes = ArrayUtils.add(abnormalNodes, i);
			}
		}
		int index = new Random().nextInt(abnormalNodes.length);
		double[] x0 = new double[graphSize];
		for (int i = 0; i < x0.length; i++) {
			x0[i] = 0.0D;
		}
		x0[abnormalNodes[index]] = 1.0D;
		return new ArrayRealVector(x0);
	}

	/**
	 * get a support of a vector
	 *
	 * @param x
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

	/**
	 * normalize a vector. if x[i] < 0.0, then x[i]:= 0.0 if x[i] > 1.0, then
	 * x[i]:= 1.0
	 * 
	 * @param x
	 *            input vector
	 * @return
	 */
	private ArrayRealVector normalize(ArrayRealVector x) {
		ArrayRealVector result = new ArrayRealVector(x);
		for (int i = 0; i < x.getDimension(); i++) {
			if (x.getEntry(i) < 0.0D) {
				result.setEntry(i, 0.0D);
			}
			if (x.getEntry(i) > 1.0D) {
				result.setEntry(i, 1.0D);
			}
		}
		return result;
	}

}
