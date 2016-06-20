package edu.albany.cs.tailApprox;

import edu.albany.cs.base.DisjointSet;
import edu.albany.cs.fastPCST.FastPCST;
import edu.albany.cs.headApprox.Digraph;
import edu.albany.cs.headApprox.DirectedEulerianCycle;
import org.apache.commons.lang3.ArrayUtils;

import java.util.*;

/**
 * Tail Approximation
 * 
 * Aglorithm 1 of
 * "http://people.csail.mit.edu/ludwigs/papers/icml15_graphsparsity.pdf" Title :
 * "A Nearly-Linear Time Framework for Graph-Structured Sparsity" Authors :
 * Ludwig Schmidt, Chinmay Hegde, and Piotr Indyk Coded by Baojian
 *
 * @author Baojian
 */
public class TailApprox {

	/** Graph G is split into 3 parts edges, edgeCosts, and prizesPi. */
	/** represents graph G */
	private final ArrayList<Integer[]> edges;
	/** edge costs c(e) (in algorithm 1 c(e) = w(e) + B / s). */
	private final ArrayList<Double> edgeCostsc;
	/** node prizes pi */
	private final ArrayList<Double> prizesPi;
	private final double C;
	private final double nu;
	/** delta */
	private final double delta;
	/** g (number of active clusters in PCST's forest) */
	private final int g; // 7.

	private int[] trueSubGraph;

	private int verboseLevel = 0;

	public F bestForest;
	public boolean valid;

	/**
	 * This is general cases (for any situation)
	 *
	 * @param edges
	 *            // Graph G
	 * @param edgeCostsc
	 *            // edge Cost c
	 * @param prizesPi
	 *            // node prizes Pi
	 * @param g
	 *            // number of connected components g
	 * @param costBudgetC
	 *            // cost Budget C
	 * @param nu
	 *            // constant parameter nu
	 * @param delta
	 *            // constant parameter delta
	 */
	public TailApprox(ArrayList<Integer[]> edges, ArrayList<Double> edgeCostsc, ArrayList<Double> prizesPi, int g,
			double costBudgetC, double nu, double delta) {

		this.edges = edges;
		this.edgeCostsc = edgeCostsc;
		this.prizesPi = prizesPi;
		this.g = g;
		this.C = costBudgetC;
		this.nu = nu;
		this.delta = delta;

		this.bestForest = run();

		/**
		 * check equation 8 : the parameter of equation 8 is cT > 1 ; cT is
		 * arbitrary and fixed constant.
		 */
		double cT = Math.sqrt(1.0D + 3.0D / (this.nu - 2.0D)); //
		double[] b = new double[prizesPi.size()];
		for (int i = 0; i < prizesPi.size(); i++) {
			b[i] = prizesPi.get(i);
		}
		double[] bS = new double[b.length];
		double[] bSPrime = new double[b.length];
		for (int i = 0; i < b.length; i++) {
			if (bestForest.nodesInF.contains(i)) {
				bS[i] = b[i];
			} else {
				bS[i] = 0.0D;
			}
			if (ArrayUtils.contains(this.trueSubGraph, i)) {
				bSPrime[i] = b[i];
			} else {
				bSPrime[i] = 0.0D;
			}
		}

		if (this.checkEqu8Valid(cT, b, bS, bSPrime)) {
			this.valid = true;
			System.out.println("result is valid");
		} else {
			this.valid = false;
			System.out.println("result of Tail approximation is invalid");
			System.exit(0);
		}
	}

	/**
	 * Note: This is for algorithm 1 only
	 *
	 * @param edges
	 * @param edgeCostsW
	 * @param z
	 * @param s
	 * @param g
	 * @param B
	 */
	public TailApprox(ArrayList<Integer[]> edges, ArrayList<Double> edgeCostsW, double[] z, int s, int g, double B,
			int[] trueSubGraph) {

		this.trueSubGraph = trueSubGraph;
		this.trueSubGraph = null; // TODO
		this.edges = edges; // 1. edges
		this.edgeCostsc = new ArrayList<Double>(); // 2. this is edge costs c
		for (double w : edgeCostsW) {
			this.edgeCostsc.add(w + B / (s * 1.0D));
		}
		this.prizesPi = new ArrayList<Double>(); // 3. prize pi

		if (z == null) {
			System.out.println("Error in PCSFTail : vector z is null ...");
			new IllegalAccessException("Error in PCSFTail : vector z is null ...");
		}
		for (int i = 0; i < z.length; i++) {
			this.prizesPi.add(z[i] * z[i]); // the prizes of pi equals to z*z in
			// algorithm 1
		}

		this.g = g;
		this.C = 2.0D * B; // C = 2.0*B in algorithm 1
		this.nu = 21; // we set nu as constant value
		this.delta = Math.min(1 / 2.0D, 1 / this.nu); // delta = min(1/2,1/v)

		this.bestForest = run();

		if (verboseLevel >= 1) {
			System.out.println("prize : " + Arrays.toString(z));
			System.out.println("prize : " + z.length);
			System.out.println("g : " + this.g + " s : " + s + " nu : " + this.nu + " " + " C : " + this.C + " delta : "
					+ this.delta);
		}

		/**
		 * try { FileWriter fileWriter = new
		 * FileWriter("./tmp/pcsfTail.txt",true) ;
		 * fileWriter.write(this.bestForest.nodesInF.size() + " "+ (2*this.nu*s
		 * + g)+"\n"); fileWriter.close() ; } catch (IOException e) {
		 * e.printStackTrace(); }
		 */

		double cT = Math.sqrt(1.0D + 3.0D / (this.nu - 2.0D)); // the parameter
																// of
		// equation 8 is cT >
		// 1 ; cT is
		// arbitrary and
		// fixed constant.
		double[] b = new double[z.length];
		for (int i = 0; i < z.length; i++) {
			b[i] = z[i];
		}
		double[] bS = new double[b.length];
		double[] bSPrime = new double[b.length];
		for (int i = 0; i < b.length; i++) {
			if (this.bestForest.nodesInF.contains(i)) {
				bS[i] = b[i];
			} else {
				bS[i] = 0.0D;
			}
			if (ArrayUtils.contains(this.trueSubGraph, i)) {
				bSPrime[i] = b[i];
			} else {
				bSPrime[i] = 0.0D;
			}
		}
		// check equation 8
		if (this.checkEqu8Valid(cT, b, bS, bSPrime)) {
			this.valid = true;
			// System.out.println("Tail approximation is valid");
		} else {
			this.valid = false;
			System.out.println("result of Tail approximation is invalid");
			System.exit(0);
		}
	}

	private F run() {
		// lambda is trade off parameter
		double minPi = this.getMin();
		double lambda0 = minPi / (2.0D * C);

		F forest = this.PCSF_GW(edges, getCostsLambda(lambda0), this.prizesPi, g);
		if (verboseLevel >= 1) {
			System.out.println("pi_min : " + minPi + " lambda0 : " + lambda0 + " c(F) : " + forest.costF + " 2*C : "
					+ 2.0D * C + " pi(\bar{F}) : " + forest.piFBar);
		}
		if ((forest.costF <= 2.0D * C) && (forest.piFBar == 0.0D || forest.piFBar <= 1e-06D)) { // check
			// this
			// further,
			// ==
			// 0.0D
			// ???
			if (forest.piFBar <= 1e-06D && forest.piFBar > 0.0D) {
				System.out.println("Error : bar{F} shoul not be nonzero ...");
				System.exit(0);
			}
			return forest;
		}
		double lambdaR = 0.0D;
		double lambdaL = 3.0D * this.getSumPrizePi();
		double epsilon = (minPi * delta) / C;
		while ((lambdaL - lambdaR) > epsilon) {
			double lambdaM = (lambdaL + lambdaR) / 2.0D;
			forest = this.PCSF_GW(edges, this.getCostsLambda(lambdaM), this.prizesPi, g);
			if ((forest.costF >= 2.0D * C) && (forest.costF <= nu * C)) {
				return forest;
			}
			if (forest.costF > nu * C) {
				lambdaR = lambdaM;
			} else {
				lambdaL = lambdaM;
			}
		}
		forest = this.PCSF_GW(edges, this.getCostsLambda(lambdaL), this.prizesPi, g);
		return forest;
	}

	/**
	 * we need to check the validation of our result. We assume the norm
	 * ||b-bSPrime|| is the minimum one. (This means you need to know which
	 * subset is optimal)
	 *
	 * @return true the equation 8 is satisfied. ; if it returns false, the
	 *         algorithm must have error(s).
	 */
	private boolean checkEqu8Valid(double cT, double[] b, double[] bS, double[] bSPrime) {

		if (this.trueSubGraph == null) { // if there is no true subgraph the
			// equation will not be checked.
			return true;
		}
		boolean flag = false;
		double leftNorm = 0.0D;
		if (b == null || bS == null || bSPrime == null || b.length == 0 || bS.length == 0 || bSPrime.length == 0) {
			return flag;
		}
		for (int i = 0; i < b.length; i++) {
			leftNorm += (b[i] - bS[i]) * (b[i] - bS[i]);
		}
		leftNorm = Math.sqrt(leftNorm);
		double rightNorm = 0.0D;
		for (int i = 0; i < b.length; i++) {
			rightNorm += (b[i] - bSPrime[i]) * (b[i] - bSPrime[i]);
		}
		rightNorm = Math.sqrt(rightNorm);
		if (leftNorm <= cT * rightNorm) {
			flag = true;
		} else {
			flag = false;
		}
		return flag;
	}

	private double getMin() {
		double minPi = Double.MAX_VALUE;
		for (int i = 0; i < this.prizesPi.size(); i++) {
			if (this.prizesPi.get(i) < minPi && this.prizesPi.get(i) > 0.0D) {
				minPi = this.prizesPi.get(i);
			}
		}
		return minPi;
	}

	private double getSumPrizePi() {
		double sumPrizePi = 0.0D;
		for (int i = 0; i < this.prizesPi.size(); i++) {
			sumPrizePi += this.prizesPi.get(i);
		}
		return sumPrizePi;
	}

	private ArrayList<Double> getCostsLambda(double lambda) {
		ArrayList<Double> cLambda = new ArrayList<Double>();
		for (int i = 0; i < this.edgeCostsc.size(); i++) {
			cLambda.add(this.edgeCostsc.get(i) * lambda);
		}
		return cLambda;
	}

	private F PCSF_GW(ArrayList<Integer[]> edges, ArrayList<Double> cLambda, ArrayList<Double> pi, int g) {

		for (double p : pi) { // check prize valid
			if (p < 0.0D) {
				new IllegalAccessException("the prize should not be negative ...");
				System.exit(0);
			}
		}

		if (!isConnected(edges)) {
			System.out.println("the graph is not connected ...");
			System.exit(0);
		}

		FastPCST pcstFast = new FastPCST(edges, pi, cLambda, FastPCST.kNoRoot, g, FastPCST.PruningMethod.kStrongPruning,
				-1);
		ArrayList<Integer> nodesInF = null;
		ArrayList<Integer> resultEdges = null;
		if (!pcstFast.run(resultEdges, nodesInF)) {
			new IllegalArgumentException("Error : Algorithm returned false. There must be an error. \n");
			System.exit(0);
		} else {
			resultEdges = pcstFast.resultEdges; // index
			nodesInF = pcstFast.resultNodes;
		}
		ArrayList<Integer[]> edgesInF = new ArrayList<Integer[]>();
		ArrayList<Double> costsInF = new ArrayList<Double>();
		if (resultEdges != null) {
			for (int i : resultEdges) {
				edgesInF.add(edges.get(i));
				costsInF.add(this.edgeCostsc.get(i));
			}
		}
		ArrayList<Double> piInF = new ArrayList<Double>();
		for (int i : nodesInF) {
			piInF.add(this.prizesPi.get(i));
		}
		double totalPrizesInG = 0.0D;
		for (int i = 0; i < this.prizesPi.size(); i++) {
			totalPrizesInG += this.prizesPi.get(i);
		}
		F forest = new F(nodesInF, edgesInF, costsInF, piInF, totalPrizesInG);
		if (forest.gamma != g) {
			System.out.println("Error : the number of trees in forest does not match with g ...");
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
	 * Forest represents the PCSF result.
	 *
	 * @author baojian
	 */
	public class F {

		public final ArrayList<Integer> nodesInF;
		public final ArrayList<Integer[]> edgesInF;
		public final ArrayList<Double> edgesCostInF;
		public final ArrayList<Double> prizePiInF;
		public final double costF;
		public final double piFBar;
		public final double prizeF;
		public final ArrayList<Tree> trees;

		public ArrayList<ArrayList<Integer[]>> pruningTrees;
		public ArrayList<Integer> pruningNodes;

		public int gamma; // The number of trees in the forest.

		public F(ArrayList<Integer> nodesInF, ArrayList<Integer[]> edgesInF, ArrayList<Double> edgesCostInF,
				ArrayList<Double> prizePiInF, double totalPrizesInG) {
			this.nodesInF = nodesInF;
			this.edgesInF = edgesInF;
			this.prizePiInF = prizePiInF;
			this.edgesCostInF = edgesCostInF;
			this.costF = getCostF();
			double result = 0.0D;
			for (double d : prizePiInF) {
				result += d;
			}
			this.prizeF = result;
			this.piFBar = totalPrizesInG - result;
			this.gamma = this.getConnectedComponents();
			this.trees = this.constructTrees();
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
			// check duplicated nodes
			HashSet<Integer> allnodes = new HashSet<Integer>();
			for (Tree tree : trees) {
				for (Integer node : tree.nodesInT) {
					if (allnodes.add(node)) { // do nothing
					} else {
						System.out.println("Error : duplicated nodes found in pruneForest ...");
						System.exit(0);
					}
				}
			}
			for (Tree tree : trees) {
				this.nodesInF.addAll(tree.nodesInT);
				this.prizePiInF.addAll(tree.prizePiInT);
				if (tree.edgesInT != null) { // check the null value
					this.edgesInF.addAll(tree.edgesInT);
					this.edgesCostInF.addAll(tree.edgesCostInT);
				}
			}
			this.costF = getCostF();
			this.prizeF = getPrizeF();
			this.piFBar = 0.0D; // be careful this will not be used.
			this.trees = constructTrees();
			this.gamma = this.getNumConnectedComponents();
			if (this.gamma != this.trees.size()) { // make sure the gamma is the
				// number of trees in this forest
				System.out.println("Error : the number of trees is not equal to gamma function ...");
				System.out.println("gamma is " + this.gamma + " is not equal to " + this.trees.size());
				System.out.println("trees size : " + trees.size());
				System.exit(0);
			}
		}

		/**
		 * @return the cost of this forest
		 */
		private double getCostF() {
			double result = 0.0D;
			for (double d : this.edgesCostInF) {
				result += d;
			}
			return result;
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
		 * @return the number of connected components in the graph
		 */
		private int getConnectedComponents() {
			DisjointSet<Integer> dis = new DisjointSet<Integer>(); // need to
																	// test
			// further TODO
			for (Integer node : nodesInF) {
				dis.makeSet(node);
			}
			for (Integer[] edge : edgesInF) {
				dis.union(edge[0], edge[1]);
			}
			return dis.numConnectedComponents;
		}

		/**
		 * construct trees for forest F
		 *
		 * @return
		 */
		public ArrayList<Tree> constructTrees() {

			DisjointSet<Integer> dis = new DisjointSet<Integer>();
			for (Integer node : nodesInF) {
				dis.makeSet(node);
			}
			for (Integer[] edge : edgesInF) {
				dis.union(edge[0], edge[1]);
			}
			HashMap<Integer, Set<Integer>> componentsMap = dis.getConnectedComponents(); // all
			// of
			// components
			// in
			// F
			ArrayList<Tree> trees = new ArrayList<Tree>();
			for (Integer key : componentsMap.keySet()) { // for each component,
															// create
				// a new tree
				Tree tree;
				Set<Integer> nodes = componentsMap.get(key); // nodes in this
																// component
				ArrayList<Integer> nodesInT = new ArrayList<Integer>(nodes); // get
																				// all
				// nodes in
				// a tree
				ArrayList<Integer[]> edgesInT = new ArrayList<Integer[]>();
				ArrayList<Double> edgesCostInT = new ArrayList<Double>();
				ArrayList<Double> prizePiInT = new ArrayList<Double>();
				for (Integer[] edge : this.edgesInF) { // get all edges in a
														// tree
					if (nodes.contains(edge[0]) || nodes.contains(edge[1])) {
						edgesInT.add(edge);
						int index = this.edgesInF.indexOf(edge);
						edgesCostInT.add(this.edgesCostInF.get(index)); // get
																		// all
																		// edge
																		// cost
						// in a tree
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

		public int[] pruningTrees(double[] y, ArrayList<Integer> tailNodes) {

			ArrayList<Integer> allNodes = new ArrayList<Integer>();
			ArrayList<ArrayList<Integer[]>> allEdges = new ArrayList<ArrayList<Integer[]>>();
			for (Tree tree : this.trees) {
				System.out.println("====================pruning tree =====================");
				HashMap<Integer, HashMap<Integer, Double>> adj = tree.adj;
				for (int i : adj.keySet()) {
					System.out.println("nei " + i + " : " + adj.get(i).keySet());
				}

				while (true) {
					boolean flag = false;
					HashMap<Integer, HashMap<Integer, Double>> updatedAdj = new HashMap<Integer, HashMap<Integer, Double>>();
					for (int k : adj.keySet()) {
						updatedAdj.put(k, new HashMap<Integer, Double>());
						for (int kk : adj.get(k).keySet()) {
							updatedAdj.get(k).put(kk, 0.0D);
						}
					}

					for (int i : adj.keySet()) {
						int nei = tree.isLeafNode(i);
						if (nei != -1) { // i is a leaf node and i is normal
											// node
							System.out.println("leaf : " + i + "y[" + i + "]=" + y[i]);
							if (y[i] == 0.0D) { // i is a leaf node and it is
												// normal node
								updatedAdj.get(nei).remove(i);
								updatedAdj.remove(i);
								flag = true;
							} else { // i is a leaf node and it is abnormal node
								int backNode = i;
								int next = nei;
								ArrayList<Integer> path = new ArrayList<Integer>();
								path.add(backNode);

								while (adj.get(next).size() == 2) {
									ArrayList<Integer> neis = new ArrayList<Integer>(adj.get(next).keySet());
									int nei0 = neis.get(0);
									int nei1 = neis.get(1);
									if (nei0 == backNode) {
										backNode = next;
										next = nei1;
									} else {
										backNode = next;
										next = nei0;
									}
									path.add(backNode);
								}
								// update adj
								boolean flag1 = false;
								for (int kk : adj.keySet()) {
									if (adj.get(kk).size() > 2) {
										flag1 = true;
										break;
									}
								}
								if (path.size() >= 3 && flag1 == true) {

									double ratio = 0.0D;
									for (int k : path) {
										if (y[k] == 0.0D) {
											ratio += 1.0D;
										}
									}
									ratio = ratio / (path.size() + 0.0D);
									System.out.println("ratio : " + ratio + " path is : " + path.toString());
									if (ratio > 0.5) {
										int lastNode = path.get(path.size() - 1);
										int nextLast = path.get(path.size() - 2);
										updatedAdj.get(lastNode).remove(nextLast);
										path.remove(path.size() - 1);
										for (int pathNode : path) {
											updatedAdj.remove(pathNode);
										}
										flag = true;
									}
								}
							}
						}
					}
					if (flag == false) {
						break;
					}
					adj = updatedAdj;
					System.out.println("end of iteration ...");
				}

				ArrayList<Integer[]> currentTree = new ArrayList<Integer[]>();
				for (int key : adj.keySet()) {
					allNodes.add(key);
					for (int k : adj.get(key).keySet()) {
						currentTree.add(new Integer[] { key, k });
					}
				}
				allEdges.add(currentTree);
			} // next tree
			int[] result = new int[allNodes.size()];
			for (int i = 0; i < allNodes.size(); i++) {
				result[i] = allNodes.get(i);
			}
			this.pruningTrees = allEdges;
			this.pruningNodes = allNodes;
			return result;
		}
	}

	/**
	 * The tree represented by edges and nodes.
	 *
	 * @author baojian
	 */
	public class Tree implements Comparable<Tree> {
		public final ArrayList<Integer> nodesInT;
		public final ArrayList<Integer[]> edgesInT; // Note : edges in T are
		// undirected. Therefore, each
		// of edge only save once.
		public final ArrayList<Double> edgesCostInT;
		public final ArrayList<Double> prizePiInT;
		public final double costTree;
		public final double prizeTree;
		public final double ratio;

		// adj for path computing
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

		public int isLeafNode(int node) {
			if (this.adj.get(node).size() == 1) {
				for (int nei : this.adj.get(node).keySet()) {
					return nei;
				}
				return -1;
			} else {
				return -1;
			}
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
}
