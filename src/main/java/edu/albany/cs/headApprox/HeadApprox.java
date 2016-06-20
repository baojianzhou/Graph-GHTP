package edu.albany.cs.headApprox;

import edu.albany.cs.base.DisjointSet;
import edu.albany.cs.fastPCST.FastPCST;
import org.apache.commons.lang3.ArrayUtils;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Head approximation
 * 
 * Aglorithm 3 of
 * "http://people.csail.mit.edu/ludwigs/papers/icml15_graphsparsity.pdf" Title :
 * "A Nearly-Linear Time Framework for Graph-Structured Sparsity" Authors :
 * Ludwig Schmidt, Chinmay Hegde, and Piotr Indyk
 *
 * @author baojian bzhou6@albany.edu
 */
public class HeadApprox {
	/** 1. Graph G */
	private ArrayList<Integer[]> edges;
	/** 2. edge costs c (corresponding to to edges) */
	private ArrayList<Double> c;
	/** 3. node prizes pi (corresponding to nodes 0,...,n-1) */
	private ArrayList<Double> pi;
	/** 4. number of connected components */
	private int g;
	/** 5. cost budget */
	private double C;
	/** 6. delta value, which is a constant value */
	private double delta;

	private double cH;
	private int[] trueSubGraph;
	private int verboseLevel = 0;

	public F bestForest;
	public boolean valid;

	public HeadApprox() {
	}

	/**
	 * A general input interface
	 * 
	 * @param edges
	 * @param c
	 * @param pi
	 * @param g
	 * @param C
	 */
	public HeadApprox(ArrayList<Integer[]> edges, ArrayList<Double> c, ArrayList<Double> pi, int g, double C) {

		this.edges = edges;
		this.c = c;
		this.pi = pi;
		this.C = C;
		this.g = g;
		/** By Theorem 11 in this paper (p.20). */
		this.delta = 1.0D / 169.0D; //
		this.cH = Math.sqrt(1.0D / 14.0D);

		this.checkInputValidation(edges, pi, c);
		this.bestForest = run();

		double[] b = new double[pi.size()];
		for (int i = 0; i < pi.size(); i++) {
			b[i] = pi.get(i);
		}

		if (this.checkEqu9Valid(b)) {
			this.valid = true;
		} else {
			System.out.println("Head approximation is invalid ...");
			System.exit(0);
			this.valid = false;
		}
	}

	/**
	 * Note : This is for algorithm 1 only
	 * 
	 * @param edges
	 * @param edgeCostsW
	 * @param b
	 * @param s
	 * @param g
	 * @param B
	 */
	public HeadApprox(ArrayList<Integer[]> edges, ArrayList<Double> edgeCostsW, double[] b, int s, int g, double B,
			int[] trueSubGraph) {
		this.edges = edges;
		/**let c(e) = w(e) + (B / s)*/
		this.c = new ArrayList<Double>();  
		for (double w : edgeCostsW) {
			this.c.add(w + B / (s * 1.0D));
		}
		/** let pi(i) = bi*bi */
		this.pi = new ArrayList<Double>(); 
		for (int i = 0; i < b.length; i++) {
			this.pi.add(b[i] * b[i]);
		}
		/** let C = 2*B */
		this.C = 2.0D * B;
		/** By Theorem 11 in this paper (p.20). */
		this.delta = 1.0D / 169.0D;
		this.cH = Math.sqrt(1.0D / 14.0D);
		this.g = g;
		this.trueSubGraph = trueSubGraph;
		/** we skip the validation of equation 9 */
		this.trueSubGraph = null;
		if (verboseLevel >= 1) {
			// ------------------------------------------------------------------------------------
			System.out.println("pi : ");
			for (int k = 0; k < this.pi.size(); k++) {
				System.out.format(",%.1f", this.pi.get(k));
			}
			System.out.println();
			System.out.println(" c(e) : " + this.c.toString());
			System.out.println("B : " + B + " C : " + this.C + " delta : " + this.delta + " cH : " + this.cH);
			// ------------------------------------------------------------------------------------
		}
		this.bestForest = run();
		if (this.checkEqu9Valid(b)) {
			this.valid = true;
		} else {
			System.out.println("Head approximation is invalid ...");
			System.exit(0);
			this.valid = false;
		}
	}

	/**
	 * check input validation, to make sure the following properties : 1. the
	 * graph is connected. 2. the graph does not have self cyclic edge
	 *
	 * @param edges
	 * @param pi
	 */
	private void checkInputValidation(ArrayList<Integer[]> edges, ArrayList<Double> pi, ArrayList<Double> edgeCosts) {
		/** check input validation */
		DisjointSet<Integer> dis = new DisjointSet<Integer>();
		HashSet<Integer> nodes = new HashSet<Integer>();
		for (Integer[] edge : edges) {
			nodes.add(edge[0]);
			nodes.add(edge[1]);
			dis.makeSet(edge[0]);
			dis.makeSet(edge[1]);
			dis.union(edge[0], edge[1]);
			if (edge[0].intValue() == edge[1].intValue()) {
				new IllegalArgumentException("bad edge : [" + edge[0] + "," + edge[1] + "]");
			}
		}
		if (dis.numConnectedComponents != 1) {
			new IllegalArgumentException("Error : the graph is not connected ...");
			System.exit(0);
		}
		if (nodes.size() != pi.size()) {
			new IllegalArgumentException("Error : edges of graph do not have whole nodes  ...");
			System.exit(0);
		}
		/** check validation of every weight of edge in the graph */
		for (double edgecost : edgeCosts) {
			if (edgecost <= 1e-9D) {
				new IllegalArgumentException("bad edge : " + edgecost);
				System.exit(0);
			}
		}
	}

	/**
	 * algorithm 3
	 *
	 * @return F a forest
	 */
	private F run() {

		double minPi = this.getMin();
		BigDecimal lambdaR = new BigDecimal((2.0D * this.C) / minPi);
		F forest = this.PCSF_GW(this.edges, this.c, this.getLambdaPi(lambdaR.doubleValue()), g);
		/** make sure we have the invariant c(Fr) > 2*C */
		if (forest.costF <= 2.0D * this.C) {
			return forest;
		}

		BigDecimal epsilon = new BigDecimal((this.delta * this.C) / (2.0D * this.getSigmaPi()));
		BigDecimal lambdaL = new BigDecimal(1.0D / (4.0D * this.getSigmaPi()));

		if (verboseLevel >= 1) {
			System.out.println("costF : " + forest.costF + " ; 2*C: " + 2.0D * this.C);
			System.out.println("min pi : " + minPi + " ; lambdaR : " + lambdaR.doubleValue());
			System.out.println("epsilon : " + epsilon.doubleValue() + " ; lambdaL : " + lambdaL.doubleValue());
		}

		int iter = 0;
		/** Binary search over the Lagrange parameter lambda */
		while (lambdaR.subtract(lambdaL).compareTo(epsilon) == 1) {
			BigDecimal lambdaM = (lambdaL.add(lambdaR)).divide(new BigDecimal(2.0D));
			forest = this.PCSF_GW(this.edges, this.c, this.getLambdaPi(lambdaM.doubleValue()), this.g);
			if (forest.costF > 2.0D * this.C) {
				lambdaR = lambdaM;
			} else {
				lambdaL = lambdaM;
			}
			if (verboseLevel >= 1) {
				iter++;
				System.out.println("iteration : " + iter);
				System.out.println("lambdaM : " + lambdaM);
				System.out.println(" epsilon : " + epsilon + "lambda R: " + lambdaR + "lambda L : " + lambdaL
						+ "lambda M : " + lambdaM);
			}
		}

		if (verboseLevel >= 1) {
			System.out
					.println("final lambdaL : " + lambdaL.doubleValue() + " final lambdaR : " + lambdaR.doubleValue());
			System.out.println("total iterations : " + iter);
		}
		F forestL = this.PCSF_GW(this.edges, this.c, this.getLambdaPi(lambdaL.doubleValue()), g);
		F forestR = this.PCSF_GW(this.edges, this.c, this.getLambdaPi(lambdaR.doubleValue()), g);
		if (verboseLevel >= 1) {
			System.out.println("# nodes in L : " + forestL.nodesInF.size() + "# edges in L :" + forestL.edgesInF.size()
					+ " # CC : " + forestL.gamma);
			System.out.println("# nodes in R : " + forestR.nodesInF.size() + "# edges in R :" + forestR.edgesInF.size()
					+ " # CC : " + forestR.gamma);
		}
		/** Prune the potentially large solution Fr */
		F forestRPrime = this.pruneForest(forestR, this.c, this.pi, this.C);
		if (forestL.prizeF >= forestRPrime.prizeF) {
			return forestL;
		} else {
			return forestRPrime;
		}
	}

	/**
	 * get total prize of this graph
	 *
	 * @return the \Sigma_{pi}
	 */
	private double getSigmaPi() {
		double result = 0.0D;
		for (double p : this.pi) {
			result += p;
		}
		return result;
	}

	/**
	 * we need to check the validation of our result. We assume the norm
	 * ||b-bSPrime|| is the minimum one. (This means you need to know which
	 * subset is optimal)
	 *
	 * @return true the equation 9 is satisfied. ; if it returns false, the
	 *         algorithm must have error(s).
	 */
	private boolean checkEqu9Valid(double[] b) {
		/** do not need to check equation 9 */
		if (this.trueSubGraph == null) {
			return true;
		}
		double[] bS = new double[b.length];
		double[] bSPrime = new double[b.length];
		for (int i = 0; i < b.length; i++) {
			if (ArrayUtils.contains(trueSubGraph, i)) {
				bSPrime[i] = b[i];
			} else {
				bSPrime[i] = 0.0D;
			}
			if (this.bestForest.nodesInF.contains(i)) {
				bS[i] = b[i];
			} else {
				bS[i] = 0.0D;
			}
		}

		boolean flag = false;
		if (bS == null || bSPrime == null || bS.length == 0 || bSPrime.length == 0) {
			return flag;
		}
		double leftNorm = 0.0D;
		for (int i = 0; i < bS.length; i++) {
			leftNorm += bS[i] * bS[i];
		}
		leftNorm = Math.sqrt(leftNorm);
		double rightNorm = 0.0D;
		for (int i = 0; i < bSPrime.length; i++) {
			rightNorm += bSPrime[i] * bSPrime[i];
		}
		rightNorm = Math.sqrt(rightNorm);
		if (leftNorm >= this.cH * rightNorm) {
			flag = true;
		} else {
			flag = false;
		}
		return flag;
	}

	/**
	 * algorithm 4 Head Approximation for the WGM : subroutine PruneForest
	 *
	 * @param forest
	 * @param edgeCosts
	 * @param prizes
	 * @param costBudget
	 * @return
	 */
	private F pruneForest(F forest, ArrayList<Double> edgeCosts, ArrayList<Double> prizes, double costBudget) {
		ArrayList<Tree> trees = forest.getDescendingSortedTrees();
		double costBudgetR = costBudget;
		ArrayList<Tree> treePrimes = new ArrayList<Tree>();
		for (int i = 0; i < forest.size(); i++) {
			Tree treePrimei;
			Tree treei = trees.get(i);
			if (costBudgetR >= treei.costTree) {
				treePrimei = treei;
				treePrimes.add(treePrimei);
				/** cost budget C^i = c(T_i) */
				costBudgetR = costBudgetR - treei.costTree;
			} else if (costBudgetR > 0) {
				treePrimei = this.pruneTree(treei, edgeCosts, prizes, costBudgetR);
				treePrimes.add(treePrimei);
				/** cost budget C^i = Cr */
				costBudgetR = 0.0D;
			} else {
				treePrimei = treei.getOneNodeMaxTree();
				/** cost budget C^i = 0 */
				treePrimes.add(treePrimei);
			}
		}
		return new F(treePrimes);
	}

	private Tree pruneTree(Tree tree, ArrayList<Double> edgeCosts, ArrayList<Double> pi, double CPrime) {

		/** T = (Vt,Et), the method of getTour has been tested */
		ArrayList<Integer> tourL = tree.getEulerTour();
		ArrayList<Double> piPrime = new ArrayList<Double>();
		HashSet<Integer> hashSet = new HashSet<Integer>();
		for (int j = 0; j < tourL.size(); j++) {
			int nodeI = tourL.get(j);
			int indexNodeI = tree.nodesInT.indexOf(nodeI);
			/** if position j is the first appearance of vj in L */
			if (hashSet.add(nodeI)) {
				if (pi.get(nodeI) != tree.prizePiInT.get(indexNodeI)) {
					System.out.println("the prize of this node is inconsistent ...");
					System.out.println(pi.get(nodeI) + "is not equal to " + tree.prizePiInT.get(indexNodeI));
					System.exit(0);
				}
				piPrime.add(pi.get(nodeI));
			} else {
				piPrime.add(0.0D);
			}
		}
		double phi = tree.prizeTree / tree.costTree;
		for (int i = 0; i < tree.nodesInT.size(); i++) {
			int nodeI = tree.nodesInT.get(i);
			if (pi.get(nodeI) != tree.prizePiInT.get(i)) {
				System.out.println("the prize of this node is inconsistent ...");
				System.out.println(pi.get(i) + "is not equal to " + tree.prizePiInT.get(tree.nodesInT.indexOf(i)));
				System.exit(0);
			}
			if (pi.get(nodeI) >= (CPrime * phi / 6.0D)) {
				return new Tree(nodeI, pi.get(i));
			}
		}

		int l = 1;
		ArrayList<Integer> Pl = new ArrayList<Integer>();
		ArrayList<Integer> PlPlus1 = new ArrayList<Integer>();
		for (int i = 0; i < tourL.size(); i++) {
			Pl.add(i);
			double CPrimeInPl = 0.0D;
			double piPrimeInPl = 0.0D;
			/** Pl has at least two nodes. It has at least one edge */
			if (Pl.size() >= 2) {
				for (int ii = 0; ii < Pl.size() - 1; ii++) {
					int Pi = tourL.get(Pl.get(ii));
					int PiPlus1 = tourL.get(Pl.get(ii + 1));
					CPrimeInPl += tree.adj.get(Pi).get(PiPlus1);
				}
			}
			/** Pl has at least one node, any path has at least one node */
			if (Pl.size() >= 1) {
				for (int ii = 0; ii < Pl.size(); ii++) {
					piPrimeInPl += piPrime.get(Pl.get(ii));
				}
			}
			if (CPrimeInPl > CPrime) {
				l++;
				PlPlus1 = Pl;
				Pl = new ArrayList<Integer>();
			} else if (piPrimeInPl >= (CPrime * phi / 6.0D)) {
				/** Pl has only one node */
				if (Pl.size() == 1) {
					int currentNode = tourL.get(Pl.get(0));
					int index = tree.nodesInT.indexOf(currentNode);
					if (pi.get(currentNode) != tree.prizePiInT.get(index)) {
						System.out.println("the prize of this node is inconsistent ...");
						System.out.println(
								pi.get(i) + "is not equal to " + tree.prizePiInT.get(tree.nodesInT.indexOf(i)));
						System.exit(0);
					}
					/** create single node tree (a tree without any edge) */
					return new Tree(currentNode, pi.get(currentNode));
				} else if (Pl.size() >= 2) {
					ArrayList<Integer> nodesInT = new ArrayList<Integer>();
					ArrayList<Integer[]> edgesInT = new ArrayList<Integer[]>();
					ArrayList<Double> edgesCostInT = new ArrayList<Double>();
					ArrayList<Double> prizePiInT = new ArrayList<Double>();
					for (int ii = 0; ii < Pl.size() - 1; ii++) {
						int Pi = Pl.get(ii);
						int PiPlus1 = Pl.get(ii + 1);
						boolean flag0 = false;
						boolean flag1 = false;
						int nodeI = tourL.get(Pi);
						int nodePlusI = tourL.get(PiPlus1);
						if (!nodesInT.contains(nodeI)) {
							flag0 = true;
							nodesInT.add(nodeI);
							prizePiInT.add(pi.get(nodeI));
							/** check pi validation */
							if (pi.get(nodeI) != tree.prizePiInT.get(tree.nodesInT.indexOf(nodeI))) {
								System.out.println("the prize of this node is inconsistent ...");
								System.out.println(
										pi.get(i) + "is not equal to " + tree.prizePiInT.get(tree.nodesInT.indexOf(i)));
								System.exit(0);
							}
						}
						if (!nodesInT.contains(nodePlusI)) {
							flag1 = true;
							nodesInT.add(nodePlusI);
							prizePiInT.add(pi.get(nodePlusI));
							/** check pi validation */
							if (pi.get(nodePlusI) != tree.prizePiInT.get(tree.nodesInT.indexOf(nodePlusI))) {
								System.out.println("the prize of this node is inconsistent ...");
								System.out.println(
										pi.get(i) + "is not equal to " + tree.prizePiInT.get(tree.nodesInT.indexOf(i)));
								System.exit(0);
							}
						}
						/** add new edge if not true */
						if ((flag0 == false) && (flag1 == false)) {
							/** do nothing as the edge already exists */
						} else {
							Integer[] edge = new Integer[] { nodeI, nodePlusI };
							double cost = tree.adj.get(nodeI).get(nodePlusI);
							edgesInT.add(edge);
							edgesCostInT.add(cost);
						}
					}
					/** return the subtree of tree T on the nodes in Pl. */
					return new Tree(nodesInT, edgesInT, edgesCostInT, prizePiInT);
				} else {
					System.out.println("Error : path Pl has at leat one node ...");
					System.exit(0);
				}
			}
		} // end for
		/** algorithm will never reach this point */
		this.mergePNodes(Pl, PlPlus1, l);
		return null;
	}

	/**
	 * This method could not be called.
	 *
	 * @param P1
	 * @param P2
	 * @param l
	 */
	private void mergePNodes(ArrayList<Integer> P1, ArrayList<Integer> P2, int l) {
		new IllegalAccessException("Error : the algorithm will never reach this point .... ");
		System.exit(0);
	}

	/**
	 * get minimum value of pi vector
	 *
	 * @return
	 */
	private double getMin() {
		double minPi = Double.MAX_VALUE;
		if (this.pi == null || this.pi.isEmpty()) {
			System.out.println("Error : pi is null ...");
			System.exit(0);
		}
		for (int i = 0; i < this.pi.size(); i++) {
			if (this.pi.get(i) < minPi && this.pi.get(i) > 0.0D) {
				minPi = this.pi.get(i);
			}
		}
		return minPi;
	}

	/**
	 * get lambda*pi vector
	 *
	 * @param lambda
	 * @return pi*lambda
	 */
	private ArrayList<Double> getLambdaPi(double lambda) {
		ArrayList<Double> piLambda = new ArrayList<Double>();
		for (int i = 0; i < this.pi.size(); i++) {
			piLambda.add(lambda * this.pi.get(i));
		}
		return piLambda;
	}

	/**
	 * The PCSF algorithm of GW pruning version
	 *
	 * @param edges
	 *            the graph G
	 * @param edgeCosts
	 *            the edge costs in graph
	 * @param pi
	 *            prizes
	 * @param g
	 *            number of connected components
	 * @return the forest of pcsf algorithm
	 */
	public F PCSF_GW(ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, ArrayList<Double> pi, int g) {

		/** check prizes are valid */
		for (double p : pi) {
			if (p < 0.0D) {
				new IllegalAccessException("the prize should not be negative ...");
				System.exit(0);
			}
		}
		/** check edgeCosts are valid */
		for (double edge : edgeCosts) {
			if (edge <= 0.0D) {
				new IllegalAccessException("the edge cost should not be non-positive ...");
				System.exit(0);
			}
		}
		/** make sure the graph is connected */
		if (!isConnected(edges)) {
			System.out.println("the graph is not connected ...");
			System.exit(0);
		}
		FastPCST pcstFast = new FastPCST(edges, pi, edgeCosts, FastPCST.kNoRoot, g,
				FastPCST.PruningMethod.kStrongPruning, -1);
		ArrayList<Integer> nodesInF = null;
		ArrayList<Integer> resultEdges = null;
		if (!pcstFast.run(resultEdges, nodesInF)) {
			new IllegalArgumentException("Error : Algorithm returned false. There must be an error. \n");
			System.exit(0);
		} else {
			/** indexes of edges */
			resultEdges = pcstFast.resultEdges;
			nodesInF = pcstFast.resultNodes;
		}
		ArrayList<Integer[]> edgesInF = new ArrayList<Integer[]>();
		ArrayList<Double> costsInF = new ArrayList<Double>();
		if (resultEdges != null) {
			for (int i : resultEdges) {
				edgesInF.add(edges.get(i));
				costsInF.add(edgeCosts.get(i));
			}
		}

		/**
		 * Note: the prizes should be the original prize not local local
		 * variable pi
		 */
		ArrayList<Double> piInF = new ArrayList<Double>();
		for (int i : nodesInF) {
			piInF.add(this.pi.get(i));
		}
		double totalPrizes = 0.0D;
		for (int i = 0; i < this.pi.size(); i++) {
			totalPrizes += this.pi.get(i);
		}

		if (verboseLevel >= 2) {
			System.out.println("number of nodes in F : " + nodesInF.size());
			System.out.println("nodesInF : " + nodesInF.toString());
			System.out.print("edges in F : ");
			for (Integer[] edge : edgesInF) {
				System.out.println("," + "[" + edge[0] + "," + edge[1] + "]");
			}
			System.out.println();
		}

		/** the nodes, prizes, edges and costs must be consistent */
		F forest = new F(nodesInF, edgesInF, costsInF, piInF, totalPrizes);
		if (verboseLevel >= 2) {
			System.out.println("forest.gamma : " + forest.gamma + " ; g : " + g);
		}
		if (forest.gamma != g) {
			System.out.println("Head Approximation Error : the number of trees in forest does not match with g ...");
			System.exit(0);
		}
		return forest;
	}

	/**
	 * Use disjoint set to check whether this graph is connected.
	 *
	 * @param edges
	 * @return true is the graph is connected , false : the graph is not
	 *         connected
	 */
	private boolean isConnected(ArrayList<Integer[]> edges) {
		DisjointSet<Integer> dis = new DisjointSet<Integer>();
		Set<Integer> nodes = new HashSet<Integer>();

		for (Integer[] edge : edges) {
			nodes.add(edge[0]);
			nodes.add(edge[1]);
		}

		for (Integer node : nodes) {
			dis.makeSet(node);
		}

		for (Integer[] edge : edges) {
			dis.union(edge[0], edge[1]);
		}

		if (dis.numConnectedComponents == 1) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * The forest represents the subset of edges found by algorithm, The forest
	 * can be seen as a combination of several trees
	 *
	 * @author baojian
	 */
	public class F {
		public final ArrayList<Integer> nodesInF;
		public final ArrayList<Integer[]> edgesInF;
		public final ArrayList<Double> edgesCostInF;
		public final ArrayList<Double> prizePiInF;
		public final double costF;
		public final double prizeF;
		public final double piFBar;
		public final int gamma;
		/** we construct corresponding trees of forest */
		public final ArrayList<Tree> trees;

		/**
		 * make sure the nodes and prizes, edges and costs are correspondingly
		 * consistent, which means the edges and nodes do not start from 0 ...
		 *
		 * @param nodesInF
		 * @param edgesInF
		 * @param edgesCostInF
		 * @param prizePiInF
		 * @param totalPrizes
		 */
		public F(ArrayList<Integer> nodesInF, ArrayList<Integer[]> edgesInF, ArrayList<Double> edgesCostInF,
				ArrayList<Double> prizePiInF, double totalPrizes) {
			this.nodesInF = nodesInF;
			this.edgesInF = edgesInF;
			this.edgesCostInF = edgesCostInF;
			this.prizePiInF = prizePiInF;
			this.costF = getCostF();
			this.prizeF = getPrizeF();
			this.trees = constructTrees();
			this.piFBar = totalPrizes - this.prizeF;
			this.gamma = this.getNumConnectedComponents();
			/** make sure gamma is the number of trees in this forest */
			if (this.gamma != this.trees.size()) {
				System.out.println("Error : the number of trees is not equal ...");
				System.exit(0);
			}
		}

		/**
		 * construct a forest from trees
		 *
		 * @param trees
		 *            the trees are combined to be a new forest
		 */
		public F(ArrayList<Tree> trees) {
			this.nodesInF = new ArrayList<Integer>();
			this.edgesInF = new ArrayList<Integer[]>();
			this.edgesCostInF = new ArrayList<Double>();
			this.prizePiInF = new ArrayList<Double>();

			if (trees == null) {
				new IllegalArgumentException("Input trees are null ...");
				System.exit(0);
			}
			/** check duplicated nodes */
			HashSet<Integer> allnodes = new HashSet<Integer>();
			for (Tree tree : trees) {
				for (Integer node : tree.nodesInT) {
					if (!allnodes.add(node)) {
						System.out.println("Error : duplicated nodes found in pruneForest ...");
						System.exit(0);
					}
				}
			}
			for (Tree tree : trees) {
				this.nodesInF.addAll(tree.nodesInT);
				this.prizePiInF.addAll(tree.prizePiInT);
				/** check the null value */
				if (tree.edgesInT != null) {
					this.edgesInF.addAll(tree.edgesInT);
					this.edgesCostInF.addAll(tree.edgesCostInT);
				}
			}
			this.costF = getCostF();
			this.prizeF = getPrizeF();
			this.piFBar = 0.0D; // be careful this will not be used.
			this.trees = constructTrees();
			this.gamma = this.getNumConnectedComponents();
			/** check gamma is the number of trees in this forest */
			if (this.gamma != this.trees.size()) {
				System.out.println("Error : the number of trees is not equal to gamma function ...");
				System.out.println("gamma is " + this.gamma + " is not equal to " + this.trees.size());
				System.out.println("trees size : " + trees.size());
				System.exit(0);
			}
		}

		/**
		 * construct trees for forest F
		 *
		 * @return
		 */
		private ArrayList<Tree> constructTrees() {

			DisjointSet<Integer> dis = new DisjointSet<Integer>();
			for (Integer node : nodesInF) {
				dis.makeSet(node);
			}
			for (Integer[] edge : edgesInF) {
				dis.union(edge[0], edge[1]);
			}
			/** all of components in forest F */
			HashMap<Integer, Set<Integer>> componentsMap = dis.getConnectedComponents();
			ArrayList<Tree> trees = new ArrayList<Tree>();
			/** for each component create a new tree */
			for (Integer key : componentsMap.keySet()) {
				Tree tree;
				/** nodes in this component */
				Set<Integer> nodes = componentsMap.get(key);
				/** all nodes in a tree */
				ArrayList<Integer> nodesInT = new ArrayList<Integer>(nodes);
				ArrayList<Integer[]> edgesInT = new ArrayList<Integer[]>();
				ArrayList<Double> edgesCostInT = new ArrayList<Double>();
				ArrayList<Double> prizePiInT = new ArrayList<Double>();
				/** all edges in a tree */
				for (Integer[] edge : this.edgesInF) {
					if (nodes.contains(edge[0]) || nodes.contains(edge[1])) {
						edgesInT.add(edge);
						int index = this.edgesInF.indexOf(edge);
						/** all edge cost in a tree */
						edgesCostInT.add(this.edgesCostInF.get(index));
					}
				}
				for (Integer node : nodesInT) {
					int index = this.nodesInF.indexOf(node);
					prizePiInT.add(this.prizePiInF.get(index));
				}
				tree = new Tree(nodesInT, edgesInT, edgesCostInT, prizePiInT);
				trees.add(tree);
			}
			return trees;
		}

		/**
		 * using disjoint set data structure to find number of connected
		 * component in a forest
		 *
		 * @return number of connected components in F
		 */
		private int getNumConnectedComponents() {
			DisjointSet<Integer> dis = new DisjointSet<Integer>();
			for (Integer node : nodesInF) {
				dis.makeSet(node);
			}
			for (Integer[] edge : edgesInF) {
				dis.union(edge[0], edge[1]);
			}
			return dis.numConnectedComponents;
		}

		/**
		 * sort the trees in this forest by descending ratio r1 > r2 > r3 > ...
		 *
		 * @return the sorted trees
		 */
		public ArrayList<Tree> getDescendingSortedTrees() {
			/** Note : with descending order */
			Collections.sort(this.trees, new Comparator<Tree>() {
				public int compare(Tree o1, Tree o2) {
					return o2.compareTo(o1);
				}
			});

			/** check duplicated nodes */
			HashSet<Integer> allnodes = new HashSet<Integer>();
			for (Tree tree : trees) {
				for (Integer node : tree.nodesInT) {
					/** do nothing */
					if (allnodes.add(node)) {
					} else {
						System.out.println("Error : duplicated nodes found in pruneForest ...");
						System.exit(0);
					}
				}
			}

			return this.trees;
		}

		/**
		 * gamma function in that paper
		 * 
		 * @return the size of the trees in Forest
		 */
		public int size() {
			return this.trees.size();
		}

		/**
		 * @return total cost of this forest
		 */
		private double getCostF() {
			double result = 0.0D;
			if (this.edgesCostInF != null) {
				for (double d : this.edgesCostInF) {
					result += d;
				}
			}
			return result;
		}

		/**
		 * @return total prizes of this forest
		 */
		private double getPrizeF() {
			double result = 0.0D;
			if (prizePiInF != null) {
				for (double d : prizePiInF) {
					result += d;
				}
			}
			return result;
		}
	}// class F

	/**
	 * The tree represented by edges and nodes.
	 *
	 * @author baojian
	 */
	public class Tree implements Comparable<Tree> {
		public final ArrayList<Integer> nodesInT;
		/** Note: edges in T are undirected. So, each of edge only save once. */
		public final ArrayList<Integer[]> edgesInT;
		public final ArrayList<Double> edgesCostInT;
		public final ArrayList<Double> prizePiInT;
		public final double costTree;
		public final double prizeTree;
		public final double ratio;

		/** adj for path computing */
		public final HashMap<Integer, HashMap<Integer, Double>> adj;

		/**
		 * Construct a new tree with more than one node
		 *
		 * @param nodesInT
		 * @param edges
		 * @param edgesCost
		 * @param prizePi
		 */
		public Tree(ArrayList<Integer> nodesInT, ArrayList<Integer[]> edges, ArrayList<Double> edgesCost,
				ArrayList<Double> prizePi) {
			this.nodesInT = nodesInT;
			this.edgesInT = edges;
			this.edgesCostInT = edgesCost;
			this.prizePiInT = prizePi;
			this.costTree = getCostTree();
			this.prizeTree = getPrizeTree();
			this.ratio = this.prizeTree / this.costTree;
			this.adj = contructAdj();
		}

		/**
		 * construct a new tree with only one node
		 *
		 * @param i
		 * @param pi
		 */
		public Tree(int i, double pi) {
			ArrayList<Integer> nodes = new ArrayList<Integer>();
			nodes.add(i);
			ArrayList<Integer[]> edges = new ArrayList<Integer[]>();
			ArrayList<Double> edgesCost = new ArrayList<Double>();
			ArrayList<Double> prizePi = new ArrayList<Double>();
			prizePi.add(pi);
			this.nodesInT = nodes;
			this.edgesInT = edges;
			this.edgesCostInT = edgesCost;
			this.prizePiInT = prizePi;
			this.costTree = getCostTree();
			this.prizeTree = getPrizeTree();
			this.ratio = this.prizeTree / this.costTree;
			this.adj = contructAdj();
		}

		/**
		 * @return adjacency list of this tree <node_i,nei(node_i)>
		 */
		private HashMap<Integer, HashMap<Integer, Double>> contructAdj() {
			HashMap<Integer, HashMap<Integer, Double>> adj = new HashMap<Integer, HashMap<Integer, Double>>();
			for (int node : this.nodesInT) {
				HashMap<Integer, Double> nei = new HashMap<Integer, Double>();
				adj.put(node, nei);
			}
			if (this.edgesCostInT != null) {
				for (int i = 0; i < this.edgesCostInT.size(); i++) {
					Integer[] edge = this.edgesInT.get(i);
					adj.get(edge[0]).put(edge[1], this.edgesCostInT.get(i));
					adj.get(edge[1]).put(edge[0], this.edgesCostInT.get(i));
				}
			}
			return adj;
		}

		/**
		 * @return the total cost of this tree
		 */
		private double getCostTree() {
			double result = 0.0D;
			if (this.edgesCostInT != null) {
				for (double d : this.edgesCostInT) {
					result += d;
				}
			}
			return result;
		}

		/**
		 * @return the total prize of this tree
		 */
		private double getPrizeTree() {
			double result = 0.0D;
			if (prizePiInT != null) {
				for (double d : prizePiInT) {
					result += d;
				}
			}
			return result;
		}

		/**
		 * @return a one node tree with this tree's maximum prize
		 */
		public Tree getOneNodeMaxTree() {
			double maximumPi = -Double.MAX_VALUE;
			int nodeID = -1;
			for (int i = 0; i < prizePiInT.size(); i++) {
				double d = prizePiInT.get(i);
				if (d > maximumPi) {
					maximumPi = d;
					nodeID = nodesInT.get(i);
				}
			}
			ArrayList<Integer> nodeInT = new ArrayList<Integer>();
			nodeInT.add(nodeID);
			ArrayList<Double> prizePi = new ArrayList<Double>();
			prizePi.add(maximumPi);
			return new Tree(nodeInT, null, null, prizePi);
		}

		/**
		 * we simulate a tree as directed graph, and each edge in the tree are
		 * simulated as two edges with different edge direction. For example, an
		 * edge [i,j] in this tree will be simulated as [i->j] and [j->i].
		 * Therefore, we can construct a connected Directed Graph. After that,
		 * we can find a tour trough depth first search algorithm.
		 *
		 * @return a Euler tour with (2|T| -1) nodes.
		 */
		public ArrayList<Integer> getEulerTour() {
			ArrayList<Integer> tour = new ArrayList<Integer>();
			Digraph digraph = new Digraph(nodesInT.size());
			ArrayList<Integer[]> directedEdges = new ArrayList<Integer[]>();
			for (Integer[] edge : edgesInT) {
				int index0 = nodesInT.indexOf(edge[0]);
				int index1 = nodesInT.indexOf(edge[1]);
				Integer[] edge1 = new Integer[] { index0, index1 };
				Integer[] edge2 = new Integer[] { index1, index0 };
				directedEdges.add(edge1);
				directedEdges.add(edge2);
			}
			for (Integer[] directedEdge : directedEdges) {
				digraph.addEdge(directedEdge[0], directedEdge[1]);
			}
			DirectedEulerianCycle di = new DirectedEulerianCycle(digraph);
			for (Iterator<Integer> iterator = di.cycle().iterator(); iterator.hasNext();) {
				int indexNode = iterator.next();
				tour.add(nodesInT.get(indexNode));
			}
			if (tour.size() != (2 * this.nodesInT.size() - 1)) {
				System.out.println("Error : the length of the tour should be (2*|V_T| -1 )");
				System.exit(0);
			}
			return tour;
		}

		// @Override
		public int compareTo(Tree o) {
			if ((this.ratio - o.ratio) < 0) {
				return -1;
			} else if ((this.ratio - o.ratio) == 0.0D) {
				return 0;
			} else {
				return 1;
			}
		}

	}// Tree class

	/**
	 * test for gamma and disjoint set
	 */
	public void testGammaDisjointSet() {
		ArrayList<Integer> nodes = new ArrayList<Integer>();
		nodes.add(1);
		nodes.add(2);
		nodes.add(3);
		nodes.add(4);
		nodes.add(5);
		ArrayList<Integer[]> edges = new ArrayList<Integer[]>();
		edges.add(new Integer[] { 1, 2 });
		edges.add(new Integer[] { 1, 3 });
		edges.add(new Integer[] { 4, 5 });
		ArrayList<Double> edgeCosts = new ArrayList<Double>();
		ArrayList<Double> prizePi = new ArrayList<Double>();
		edgeCosts.add(1.0D);
		edgeCosts.add(2.0D);
		edgeCosts.add(3.0D);
		prizePi.add(1.1D);
		prizePi.add(1.2D);
		prizePi.add(1.3D);
		prizePi.add(1.4D);
		prizePi.add(1.5D);
		F f = new F(nodes, edges, edgeCosts, prizePi, 0.0D);
		System.out.println(f.gamma);
	}

	/**
	 * test the tour algorihtm
	 */
	public void testTOUR() {
		ArrayList<Integer> nodesInT = new ArrayList<Integer>();
		nodesInT.add(12);
		nodesInT.add(13);
		nodesInT.add(14);
		nodesInT.add(21);
		nodesInT.add(129);
		nodesInT.add(30);
		ArrayList<Integer[]> edges = new ArrayList<Integer[]>();
		edges.add(new Integer[] { 12, 13 });
		edges.add(new Integer[] { 12, 14 });
		edges.add(new Integer[] { 30, 14 });
		edges.add(new Integer[] { 21, 14 });
		edges.add(new Integer[] { 129, 14 });
		ArrayList<Double> edgesCost = new ArrayList<Double>();
		edgesCost.add(1.0);
		edgesCost.add(2.0);
		edgesCost.add(3.0);
		edgesCost.add(4.0);
		edgesCost.add(5.0);
		ArrayList<Double> prizePi = new ArrayList<Double>();
		prizePi.add(1.0);
		prizePi.add(2.0);
		prizePi.add(3.0);
		prizePi.add(4.0);
		prizePi.add(5.0);
		prizePi.add(6.0);
		Tree tree = new Tree(nodesInT, edges, edgesCost, prizePi);
		for (int i : tree.getEulerTour()) {
			System.out.print(i + " ");
		}
		System.out.println();
	}

	public static void main(String args[]) {
		HeadApprox head = new HeadApprox();
		head.testTOUR();
		head.testGammaDisjointSet();
	}
}
