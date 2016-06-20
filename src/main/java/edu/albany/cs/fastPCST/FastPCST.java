package edu.albany.cs.fastPCST;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * FastPCST algorithm based on the following paper
 * 
 * Hegde, Chinmay, Piotr Indyk, and Ludwig Schmidt.
 * "A fast, adaptive variant of the goemans-williamson scheme for the prize-collecting steiner tree problem."
 * Workshop of the 11th DIMACS Implementation Challenge. URL https://people.
 * csail. mit. edu/ludwigs/papers/dimacs14_fastpcst. pdf. 2014.
 * 
 * There is a C++ version https://github.com/ludwigschmidt/pcst-fast.git
 * 
 * I appreciate Ludwig Schmidt provided me the c++ code and helped me a lot when
 * I implemented this java version.
 * 
 * This java version has been tested on all of the DIMACS data. The
 * corresponding dataset can be found in https://github.com/ls-cwi/heinz.git
 * 
 *
 * @author baojian bzhou6@albany.edu
 */
public class FastPCST {

	/** defined parameters for fast pcst */

	/** public heap_buffer (maybe do not need this) */
	private ArrayList<PairingHeap.Node> pairingHeapBuffer = new ArrayList<PairingHeap.Node>();
	/** edge parts information */
	private ArrayList<EdgePart> edgeParts = new ArrayList<EdgePart>();
	private ArrayList<EdgeInfo> edgeInfo = new ArrayList<EdgeInfo>();
	/** cluster information */
	private ArrayList<Cluster> clusters = new ArrayList<Cluster>();
	private ArrayList<InactiveMergeEvent> inactiveMergeEvents = new ArrayList<InactiveMergeEvent>();
	private PriorityQueue clustersDeactivation = new PriorityQueue();
	private PriorityQueue clustersNextEdgeEvent = new PriorityQueue();
	private double currentTime;
	private double eps;

	/** marks whether a node survives simple pruning */
	private ArrayList<Boolean> nodeGood = new ArrayList<Boolean>();
	private ArrayList<Boolean> nodeDeleted = new ArrayList<Boolean>();
	private ArrayList<Integer> phase2Result = new ArrayList<Integer>();
	private ArrayList<Pair<Integer, Double>> pathCompressionVisited = new ArrayList<Pair<Integer, Double>>();
	private ArrayList<Integer> clusterQueue = new ArrayList<Integer>();
	private ArrayList<ArrayList<Pair<Integer, Double>>> phase3Neighbors = new ArrayList<ArrayList<Pair<Integer, Double>>>();

	/** variables for strong pruning */
	private ArrayList<Integer> finalComponentLabel = new ArrayList<Integer>();
	private ArrayList<ArrayList<Integer>> finalComponents = new ArrayList<ArrayList<Integer>>();
	private int rootComponentIndex;
	private ArrayList<Pair<Integer, Double>> strongPruningParent = new ArrayList<Pair<Integer, Double>>();
	private ArrayList<Double> strongPruningPayoff = new ArrayList<Double>();
	private ArrayList<Pair<Boolean, Integer>> stack = new ArrayList<Pair<Boolean, Integer>>();
	private ArrayList<Integer> stack2 = new ArrayList<Integer>();

	/** used in the input of constructor */
	public static final int kNoRoot = -1;
	private ArrayList<Integer[]> edges;
	private ArrayList<Double> costs;
	private ArrayList<Double> prizes;
	private int root;
	private int targetNumActiveClusters;
	private PruningMethod pruning;
	private int verbosity_level;

	/** results of algorithm */
	public FastPCST.Statistics stats = new FastPCST.Statistics();
	public ArrayList<Integer> resultNodes;
	public ArrayList<Integer> resultEdges;

	/** reference parameters */
	private double next_edge_time = 0.0d;
	private int next_edge_cluster_index = 0;
	private int next_edge_part_index = -1;
	private double next_cluster_time;
	private int next_cluster_index;
	private double sum_current_edge_part;
	private double current_finished_moat_sum;
	private int current_cluster_index;
	private double sum_other_edge_part;
	private int other_cluster_index;
	private double other_finished_moat_sum;
	private Cluster current_cluster;
	private Cluster other_cluster;
	private EdgePart next_edge_part;
	private EdgePart other_edge_part;
	private ArrayList<Integer> phase1_result;

	/** reference parameters */
	private ArrayList<Integer> build_phase1_node_set_second_P;
	private ArrayList<Integer> build_phase1_node_set_first_P;
	private ArrayList<Integer> build_phase2_node_set_first_P;
	private ArrayList<Integer> build_phase3_node_set_first_P;

	/**
	 * @param edges
	 * @param prizes
	 * @param costs
	 * @param root
	 * @param target_num_active_clusters
	 * @param pruningMethod
	 * @param verbostiy_level
	 */
	public FastPCST(ArrayList<Integer[]> edges, ArrayList<Double> prizes, ArrayList<Double> costs, int root,
			int target_num_active_clusters, PruningMethod pruningMethod, int verbostiy_level) {

		this.edges = edges;
		this.prizes = prizes;
		this.costs = costs;
		this.root = root;
		this.targetNumActiveClusters = target_num_active_clusters;
		this.pruning = pruningMethod;
		this.verbosity_level = verbostiy_level;

		if (!checkValidInput()) {
			System.out.println("Error : input parameters are not valid..");
			System.exit(0);
		}

		edgeParts = new ArrayList<EdgePart>();
		edgeParts = resize(edgeParts, 2 * this.edges.size(), new EdgePart());
		nodeDeleted = new ArrayList<Boolean>();
		nodeDeleted = resize(nodeDeleted, prizes.size(), false);
		edgeInfo = new ArrayList<EdgeInfo>();
		edgeInfo = resize(edgeInfo, this.edges.size(), new EdgeInfo());

		for (int ii = 0; ii < edgeInfo.size(); ++ii) {
			edgeInfo.get(ii).inactiveMergeEvent = -1;
		}
		currentTime = 0.0;
		eps = 1e-10; // set to min input value / 2.0?
		// initialize clusters clusters_deactivation
		for (int ii = 0; ii < this.prizes.size(); ++ii) {
			Cluster cluster = new Cluster(pairingHeapBuffer);
			clusters.add(cluster);
			clusters.get(ii).active = (ii != root);
			clusters.get(ii).active_start_time = 0.0;
			clusters.get(ii).active_end_time = -1.0;
			if (ii == root) {
				clusters.get(ii).active_end_time = 0.0;
			}
			clusters.get(ii).merged_into = -1;
			clusters.get(ii).prize_sum = prizes.get(ii);
			clusters.get(ii).subcluster_moat_sum = 0.0;
			clusters.get(ii).moat = 0.0;
			clusters.get(ii).contains_root = (ii == root);
			clusters.get(ii).skip_up = -1;
			clusters.get(ii).skip_up_sum = 0.0;
			clusters.get(ii).merged_along = -1;
			clusters.get(ii).child_cluster_1 = -1;
			clusters.get(ii).child_cluster_2 = -1;
			clusters.get(ii).necessary = false;
			if (clusters.get(ii).active) {
				clustersDeactivation.insert(prizes.get(ii), ii);
			}
		}
		// split edge into two parts
		for (int ii = 0; ii < this.edges.size(); ++ii) {
			int uu = this.edges.get(ii)[0];
			int vv = this.edges.get(ii)[1];
			double cost = costs.get(ii);
			EdgePart uu_part = edgeParts.get(2 * ii);
			EdgePart vv_part = edgeParts.get(2 * ii + 1);
			Cluster uu_cluster = clusters.get(uu);
			Cluster vv_cluster = clusters.get(vv);
			uu_part.deleted = false;
			vv_part.deleted = false;
			if (uu_cluster.active && vv_cluster.active) {
				double event_time = cost / 2.0d;
				uu_part.nextEventVal = event_time;
				vv_part.nextEventVal = event_time;
			} else if (uu_cluster.active) {
				uu_part.nextEventVal = cost;
				vv_part.nextEventVal = 0.0;
			} else if (vv_cluster.active) {
				uu_part.nextEventVal = 0.0;
				vv_part.nextEventVal = cost;
			} else {
				uu_part.nextEventVal = 0.0;
				vv_part.nextEventVal = 0.0;
			}
			uu_part.heapNode = uu_cluster.edge_parts.insert(uu_part.nextEventVal, 2 * ii);
			vv_part.heapNode = vv_cluster.edge_parts.insert(vv_part.nextEventVal, (2 * ii + 1));
			clusters.set(uu, uu_cluster);
			clusters.set(vv, vv_cluster);
			edgeParts.set(2 * ii, uu_part);
			edgeParts.set(2 * ii + 1, vv_part);
		}

		// initialize clusters_next_edge_event
		for (int ii = 0; ii < prizes.size(); ++ii) {
			if (clusters.get(ii).active) {
				if (!clusters.get(ii).edge_parts.is_empty()) {
					double val = 0.0;
					int edge_part = -1;
					clusters.get(ii).edge_parts.get_min(val, edge_part);
					val = clusters.get(ii).edge_parts.get_min_firstP;
					edge_part = clusters.get(ii).edge_parts.get_min_secondP;
					clustersNextEdgeEvent.insert(val, ii);
				}
			}
		}
	}// constructor

	private boolean checkValidInput() {
		for (Double prize : this.prizes) {
			if (prize < 0.0D) {
				System.out.println("Error : input parameter prizes are not valid..");
				return false;
			}
		}
		for (Double cost : this.costs) {
			if (cost < 0.0D) {
				System.out.println("Error : input parameter costs are not valid..");
				return false;
			}
		}
		return true;
	}

	private int get_other_edge_part_index(int edge_part_index) {
		if (edge_part_index % 2 == 0) {
			return (edge_part_index + 1);
		} else {
			return (edge_part_index - 1);
		}
	}

	public FastPCST.PruningMethod parse_pruning_method(String input) {
		PruningMethod result = PruningMethod.kUnknownPruning;
		String input_lower = input.toLowerCase();
		input = input_lower;
		if (input.equals("none")) {
			result = PruningMethod.kNoPruning;
		} else if (input.equals("simple")) {
			result = PruningMethod.kSimplePruning;
		} else if (input.equals("gw")) {
			result = PruningMethod.kGWPruning;
		} else if (input.equals("strong")) {
			result = PruningMethod.kStrongPruning;
		}
		return result;
	}

	public FastPCST.Statistics get_statistics() {
		return this.stats;
	}

	public void get_next_edge_event(Double next_time, Integer next_cluster_index, Integer next_edge_part_index) {
		if (clustersNextEdgeEvent.is_empty()) {
			next_time = Double.POSITIVE_INFINITY;
			next_cluster_index = -1;
			next_edge_part_index = -1;
			this.next_edge_time = next_time;
			this.next_edge_cluster_index = next_cluster_index;
			this.next_edge_part_index = next_edge_part_index;
		}
		clustersNextEdgeEvent.get_min(next_time, next_cluster_index);
		next_time = clustersNextEdgeEvent.get_min_firstP;
		next_cluster_index = clustersNextEdgeEvent.get_min_secondP;
		clusters.get(next_cluster_index).edge_parts.get_min(next_time, next_edge_part_index);
		next_time = clusters.get(next_cluster_index).edge_parts.get_min_firstP;
		next_edge_part_index = clusters.get(next_cluster_index).edge_parts.get_min_secondP;
		this.next_edge_time = next_time;
		this.next_edge_cluster_index = next_cluster_index;
		this.next_edge_part_index = next_edge_part_index;
	}

	public void remove_next_edge_event(int next_cluster_index) {
		clustersNextEdgeEvent.delete_element(next_cluster_index);
		double tmp_value = 0d;
		int tmp_edge_part = 1;
		clusters.get(next_cluster_index).edge_parts.delete_min(tmp_value, tmp_edge_part);
		tmp_value = clusters.get(next_cluster_index).edge_parts.delete_min_firstP;
		tmp_edge_part = clusters.get(next_cluster_index).edge_parts.delete_min_secondP;
		if (!clusters.get(next_cluster_index).edge_parts.is_empty()) {
			clusters.get(next_cluster_index).edge_parts.get_min(tmp_value, tmp_edge_part);
			tmp_value = clusters.get(next_cluster_index).edge_parts.get_min_firstP;
			tmp_edge_part = clusters.get(next_cluster_index).edge_parts.get_min_secondP;
			clustersNextEdgeEvent.insert(tmp_value, next_cluster_index);
		}
	}

	public void get_next_cluster_event(Double next_cluster_time, Integer next_cluster_index) {
		if (clustersDeactivation.is_empty()) {
			next_cluster_time = Double.POSITIVE_INFINITY;
			next_cluster_index = -1;
			this.next_cluster_time = next_cluster_time;
			this.next_cluster_index = next_cluster_index;
			return;
		}
		clustersDeactivation.get_min(next_cluster_time, next_cluster_index);
		next_cluster_time = clustersDeactivation.get_min_firstP;
		next_cluster_index = clustersDeactivation.get_min_secondP;
		this.next_cluster_time = clustersDeactivation.get_min_firstP;
		this.next_cluster_index = clustersDeactivation.get_min_secondP;
	}

	public void remove_next_cluster_event() {
		Double tmp_value = 0d;
		Integer tmp_cluster = 1;
		clustersDeactivation.delete_min(tmp_value, tmp_cluster);
		tmp_value = clustersDeactivation.delete_min_firstP;
		tmp_cluster = clustersDeactivation.delete_min_secondP;
	}

	public void get_sum_on_edge_part(int edge_part_index, double total_sum, Double finished_moat_sum,
			Integer cur_cluster_index, int flag) {
		int endpoint = edges.get(edge_part_index / 2)[0];
		if (edge_part_index % 2 == 1) {
			endpoint = edges.get(edge_part_index / 2)[1];
		}
		total_sum = 0.0;
		cur_cluster_index = endpoint;
		pathCompressionVisited = new ArrayList<Pair<Integer, Double>>();
		while (clusters.get(cur_cluster_index).merged_into != -1) {
			pathCompressionVisited.add(new Pair<Integer, Double>(cur_cluster_index, total_sum));
			if (clusters.get(cur_cluster_index).skip_up >= 0) {
				total_sum += clusters.get(cur_cluster_index).skip_up_sum;
				cur_cluster_index = clusters.get(cur_cluster_index).skip_up;
			} else {
				total_sum += clusters.get(cur_cluster_index).moat;
				cur_cluster_index = clusters.get(cur_cluster_index).merged_into;
			}
		}
		for (int ii = 0; ii < pathCompressionVisited.size(); ++ii) {
			int visited_cluster_index = pathCompressionVisited.get(ii).getFirst();
			double visited_sum = pathCompressionVisited.get(ii).getSecond();
			clusters.get(visited_cluster_index).skip_up = cur_cluster_index;
			clusters.get(visited_cluster_index).skip_up_sum = total_sum - visited_sum;
		}
		if (clusters.get(cur_cluster_index).active) {
			finished_moat_sum = total_sum;
			total_sum += currentTime - clusters.get(cur_cluster_index).active_start_time;
		} else {
			total_sum += clusters.get(cur_cluster_index).moat;
			finished_moat_sum = total_sum;
		}
		if (flag == 1) {
			sum_current_edge_part = total_sum;
			current_finished_moat_sum = finished_moat_sum;
			current_cluster_index = cur_cluster_index;
		} else {
			sum_other_edge_part = total_sum;
			other_finished_moat_sum = finished_moat_sum;
			other_cluster_index = cur_cluster_index;
		}
	}

	public void mark_nodes_as_good(int start_cluster_index) {
		clusterQueue = new ArrayList<Integer>();
		int queue_index = 0;
		clusterQueue.add(start_cluster_index);
		while (queue_index < (int) clusterQueue.size()) {
			int cur_cluster_index = clusterQueue.get(queue_index);
			queue_index += 1;
			if (clusters.get(cur_cluster_index).merged_along >= 0) {
				clusterQueue.add(clusters.get(cur_cluster_index).child_cluster_1);
				clusterQueue.add(clusters.get(cur_cluster_index).child_cluster_2);
			} else {
				nodeGood.set(cur_cluster_index, true);
			}
		}
	}

	public void mark_clusters_as_necessary(int start_cluster_index) {
		int cur_cluster_index = start_cluster_index;
		while (!clusters.get(cur_cluster_index).necessary) {
			clusters.get(cur_cluster_index).necessary = true;
			if (clusters.get(cur_cluster_index).merged_into >= 0) {
				cur_cluster_index = clusters.get(cur_cluster_index).merged_into;
			} else {
				return;
			}
		}
	}

	public void mark_nodes_as_deleted(int start_node_index, int parent_node_index) {
		nodeDeleted.set(start_node_index, true);
		clusterQueue = new ArrayList<Integer>();
		int queue_index = 0;
		clusterQueue.add(start_node_index);
		while (queue_index < clusterQueue.size()) {
			int cur_node_index = clusterQueue.get(queue_index);
			queue_index += 1;
			for (int ii = 0; ii < phase3Neighbors.get(cur_node_index).size(); ++ii) {
				int next_node_index = phase3Neighbors.get(cur_node_index).get(ii).getFirst();
				if (next_node_index == parent_node_index) {
					continue;
				}
				if (nodeDeleted.get(next_node_index)) {
					continue; // should never happen
				}
				nodeDeleted.set(next_node_index, true);
				clusterQueue.add(next_node_index);
			}
		}
	}

	public boolean run(ArrayList<Integer> result_nodes, ArrayList<Integer> result_edges) {

		result_nodes = new ArrayList<Integer>();
		result_edges = new ArrayList<Integer>();
		if (root >= 0 && targetNumActiveClusters > 0) {
			System.out.println("Error: target_num_active_clusters must be 0 in the rooted case.\n");
			System.exit(0);
			return false;
		}
		phase1_result = new ArrayList<Integer>();
		int num_active_clusters = prizes.size();
		if (root >= 0) {
			num_active_clusters -= 1;
		}
		// growth phase
		while (num_active_clusters > targetNumActiveClusters) {
			if (verbosity_level >= 2) {
				System.out.print("-----------------------------------------\n");
			}
			next_edge_time = 0.0;
			next_edge_cluster_index = 0;
			next_edge_part_index = -1;
			get_next_edge_event(next_edge_time, next_edge_cluster_index, next_edge_part_index);
			get_next_cluster_event(next_cluster_time, next_cluster_index);
			if (verbosity_level >= 2) {
				System.out.format("Next edge event: time %6f, cluster %d, part %d\n", next_edge_time,
						next_edge_cluster_index, next_edge_part_index);
				System.out.format("Next cluster event: time %6f, cluster %d\n", next_cluster_time, next_cluster_index);
			}
			// edge is tight
			if (next_edge_time < next_cluster_time) {
				if (verbosity_level >= 3) {
					System.out.format("next_edge_time : %6f ; next_cluster_time : %6f\n", next_edge_time,
							next_cluster_time);
				}
				stats.total_num_edge_events += 1;
				currentTime = next_edge_time;
				remove_next_edge_event(next_edge_cluster_index);
				if (edgeParts.get(next_edge_part_index).deleted) {
					stats.num_deleted_edge_events += 1;
					if (verbosity_level >= 2) {
						System.out.format("Edge part %d already deleted, nothing to do\n", next_edge_part_index);
					}
					// System.out.println("PCSF Test F:" +
					// stats.total_num_edge_events) ;
					continue;
				}
				// collect all the relevant information about the edge parts
				int other_edge_part_index = get_other_edge_part_index(next_edge_part_index);
				double current_edge_cost = costs.get(next_edge_part_index / 2);
				sum_current_edge_part = 0.0d;
				current_cluster_index = -1;
				current_finished_moat_sum = -1;
				get_sum_on_edge_part(next_edge_part_index, sum_current_edge_part, current_finished_moat_sum,
						current_cluster_index, 1);
				if (verbosity_level >= 3) {
					System.out.format(
							"next_edge_part_index : %d ; sum_current_edge_part : %6f ; current_finished_moat_sum : %6f ; current_cluster_index : %d\n",
							next_edge_part_index, sum_current_edge_part, current_finished_moat_sum,
							current_cluster_index);
				}
				sum_other_edge_part = 0.0d;
				other_cluster_index = -1;
				other_finished_moat_sum = -1d;
				get_sum_on_edge_part(other_edge_part_index, sum_other_edge_part, other_finished_moat_sum,
						other_cluster_index, 2);
				if (verbosity_level >= 3) {
					System.out.format(
							"other_edge_part_index : %d ; sum_other_edge_part : %6f ; other_finished_moat_sum : %6f ; other_cluster_index : %d\n",
							other_edge_part_index, sum_other_edge_part, other_finished_moat_sum, other_cluster_index);
				}
				double remainder = current_edge_cost - sum_current_edge_part - sum_other_edge_part;
				current_cluster = clusters.get(current_cluster_index);
				other_cluster = clusters.get(other_cluster_index);
				next_edge_part = edgeParts.get(next_edge_part_index);
				other_edge_part = edgeParts.get(other_edge_part_index);
				if (verbosity_level >= 2) {
					System.out.format(
							"Edge event at time %6f, current edge part %d (cluster %d), other edge part %d (cluster %d)\n",
							currentTime, next_edge_part_index, current_cluster_index, other_edge_part_index,
							other_cluster_index);
					System.out.format("Sum current part %6f, other part %6f, total length %6f, remainder %6f\n",
							sum_current_edge_part, sum_other_edge_part, current_edge_cost, remainder);
				}
				if (verbosity_level >= 3) {
					System.out.format("current_cluster_index : %d ; other_cluster_index : %d\n", current_cluster_index,
							other_cluster_index);
				}
				if (current_cluster_index == other_cluster_index) {
					stats.num_merged_edge_events += 1;
					if (verbosity_level >= 2) {
						System.out.println("Clusters already merged, ignoring edge\n");
					}
					edgeParts.get(other_edge_part_index).deleted = true;
					continue;
				}
				// merge two clusters
				if (remainder < eps * current_edge_cost) {
					stats.total_num_merge_events += 1;
					phase1_result.add(next_edge_part_index / 2);
					edgeParts.get(other_edge_part_index).deleted = true;
					int new_cluster_index = clusters.size();
					clusters.add(new Cluster(pairingHeapBuffer));
					Cluster new_cluster = clusters.get(new_cluster_index);
					current_cluster = clusters.get(current_cluster_index);
					other_cluster = clusters.get(other_cluster_index);
					if (verbosity_level >= 2) {
						System.out.format("Merge %d and %d into %d\n", current_cluster_index, other_cluster_index,
								new_cluster_index);
					}
					new_cluster.moat = 0.0;
					new_cluster.prize_sum = current_cluster.prize_sum + other_cluster.prize_sum;
					new_cluster.subcluster_moat_sum = current_cluster.subcluster_moat_sum
							+ other_cluster.subcluster_moat_sum;
					new_cluster.contains_root = current_cluster.contains_root || other_cluster.contains_root;
					new_cluster.active = !new_cluster.contains_root;
					new_cluster.merged_along = next_edge_part_index / 2;
					new_cluster.child_cluster_1 = current_cluster_index;
					new_cluster.child_cluster_2 = other_cluster_index;
					new_cluster.necessary = false;
					new_cluster.skip_up = -1;
					new_cluster.skip_up_sum = 0.0;
					new_cluster.merged_into = -1;
					current_cluster.active = false;
					current_cluster.active_end_time = currentTime + remainder;
					current_cluster.merged_into = new_cluster_index;
					current_cluster.moat = current_cluster.active_end_time - current_cluster.active_start_time;
					clustersDeactivation.delete_element(current_cluster_index);
					num_active_clusters -= 1;
					if (!current_cluster.edge_parts.is_empty()) {
						clustersNextEdgeEvent.delete_element(current_cluster_index);
					}
					// merge with active or inactive cluster
					if (other_cluster.active) {
						stats.num_active_active_merge_events += 1;
						other_cluster.active = false;
						other_cluster.active_end_time = currentTime + remainder;
						other_cluster.moat = other_cluster.active_end_time - other_cluster.active_start_time;
						clustersDeactivation.delete_element(other_cluster_index);
						if (!other_cluster.edge_parts.is_empty()) {
							clustersNextEdgeEvent.delete_element(other_cluster_index);
						}
						num_active_clusters -= 1;
					} else {
						stats.num_active_inactive_merge_events += 1;
						if (!other_cluster.contains_root) {
							double edge_event_update_time = currentTime + remainder - other_cluster.active_end_time;
							other_cluster.edge_parts.add_to_heap(edge_event_update_time);
							inactiveMergeEvents.add(new InactiveMergeEvent());
							InactiveMergeEvent merge_event = inactiveMergeEvents.get(inactiveMergeEvents.size() - 1);
							merge_event.active_cluster_index = current_cluster_index;
							merge_event.inactive_cluster_index = other_cluster_index;
							int active_node_part = edges.get(next_edge_part_index / 2)[0];
							int inactive_node_part = edges.get(next_edge_part_index / 2)[1];
							if ((next_edge_part_index % 2) == 1) {
								int tmp = active_node_part;
								active_node_part = inactive_node_part;
								inactive_node_part = tmp;
							}
							merge_event.active_cluster_node = active_node_part;
							merge_event.inactive_cluster_node = inactive_node_part;
							edgeInfo.get(next_edge_part_index / 2).inactiveMergeEvent = inactiveMergeEvents.size() - 1;
						}
					}
					other_cluster.merged_into = new_cluster_index;
					new_cluster.edge_parts = current_cluster.edge_parts.meld(current_cluster.edge_parts,
							other_cluster.edge_parts);
					new_cluster.subcluster_moat_sum += current_cluster.moat;
					new_cluster.subcluster_moat_sum += other_cluster.moat;
					if (new_cluster.active) {// new_cluster is inactive if new
												// cluster contains root
						new_cluster.active_start_time = currentTime + remainder;
						double becoming_inactive_time = currentTime + remainder + new_cluster.prize_sum
								- new_cluster.subcluster_moat_sum;
						clustersDeactivation.insert(becoming_inactive_time, new_cluster_index);
						if (!new_cluster.edge_parts.is_empty()) {
							double tmp_val = 0.0d;
							int tmp_index = -1;
							new_cluster.edge_parts.get_min(tmp_val, tmp_index);
							tmp_val = new_cluster.edge_parts.get_min_firstP;
							tmp_index = new_cluster.edge_parts.get_min_secondP;
							clustersNextEdgeEvent.insert(tmp_val, new_cluster_index);
						}
						num_active_clusters += 1;
					}
					// do not merge two clusters, but the edge still
					// tight.(since the other edge part is not tight)
				} else if (other_cluster.active) {
					stats.total_num_edge_growth_events += 1;
					stats.num_active_active_edge_growth_events += 1;
					double next_event_time = currentTime + remainder / 2.0;
					next_edge_part.nextEventVal = sum_current_edge_part + remainder / 2.0;
					if (!current_cluster.edge_parts.is_empty()) {
						clustersNextEdgeEvent.delete_element(current_cluster_index);
					}
					next_edge_part.heapNode = current_cluster.edge_parts.insert(next_event_time, next_edge_part_index);
					double tmp_val = -1.0;
					int tmp_index = -1;
					current_cluster.edge_parts.get_min(tmp_val, tmp_index);
					tmp_val = current_cluster.edge_parts.get_min_firstP;
					tmp_index = current_cluster.edge_parts.get_min_secondP;
					clustersNextEdgeEvent.insert(tmp_val, current_cluster_index);
					clustersNextEdgeEvent.delete_element(other_cluster_index);
					other_cluster.edge_parts.decrease_key(other_edge_part.heapNode,
							other_cluster.active_start_time + other_edge_part.nextEventVal - other_finished_moat_sum,
							next_event_time);
					other_cluster.edge_parts.get_min(tmp_val, tmp_index);
					tmp_val = other_cluster.edge_parts.get_min_firstP;
					tmp_index = other_cluster.edge_parts.get_min_secondP;
					clustersNextEdgeEvent.insert(tmp_val, other_cluster_index);
					other_edge_part.nextEventVal = sum_other_edge_part + remainder / 2.0;
					if (verbosity_level >= 2) {
						System.out.format("Added new event at time %6f\n", next_event_time);
					}
					// generate new edge events
				} else {
					stats.total_num_edge_growth_events += 1;
					stats.num_active_inactive_edge_growth_events += 1;
					double next_event_time = currentTime + remainder;
					next_edge_part.nextEventVal = current_edge_cost - other_finished_moat_sum;
					if (!current_cluster.edge_parts.is_empty()) {
						clustersNextEdgeEvent.delete_element(current_cluster_index);
					}
					next_edge_part.heapNode = current_cluster.edge_parts.insert(next_event_time, next_edge_part_index);
					double tmp_val = -1.0;
					int tmp_index = -1;
					current_cluster.edge_parts.get_min(tmp_val, tmp_index);
					tmp_val = current_cluster.edge_parts.get_min_firstP;
					tmp_index = current_cluster.edge_parts.get_min_secondP;
					clustersNextEdgeEvent.insert(tmp_val, current_cluster_index);
					other_cluster.edge_parts.decrease_key(other_edge_part.heapNode,
							other_cluster.active_end_time + other_edge_part.nextEventVal - other_finished_moat_sum,
							other_cluster.active_end_time);
					other_edge_part.nextEventVal = other_finished_moat_sum;
					if (verbosity_level >= 2) {
						System.out.format("Added new event at time %6f and event for inactive edge part\n",
								next_event_time);
					}
				}
				// cluster is tight
			} else {
				stats.num_cluster_events += 1;
				currentTime = next_cluster_time;
				remove_next_cluster_event();
				Cluster cur_cluster = clusters.get(next_cluster_index);
				cur_cluster.active = false;
				cur_cluster.active_end_time = currentTime;
				cur_cluster.moat = cur_cluster.active_end_time - cur_cluster.active_start_time;
				if (!cur_cluster.edge_parts.is_empty()) {
					clustersNextEdgeEvent.delete_element(next_cluster_index);
				}
				num_active_clusters -= 1;
				if (verbosity_level >= 2) {
					System.out.format("Cluster deactivation: cluster %d at time %6f (moat size %6f)\n",
							next_cluster_index, currentTime, cur_cluster.moat);
				}
			}
		} // while(num_active_clusters > target_num_active_clusters)

		if (verbosity_level >= 1) {
			System.out.format("Finished GW clustering: final event time %6f, number of edge events %d\n", currentTime,
					stats.total_num_edge_events);
		}
		nodeGood = new ArrayList<Boolean>();
		for (int i = 0; i < prizes.size(); i++) {
			nodeGood.add(false);
		}
		if (root >= 0) {
			System.out.println("has not checked");
			System.exit(0);
			// find the root cluster
			for (int ii = 0; ii < clusters.size(); ++ii) {
				if (clusters.get(ii).contains_root && clusters.get(ii).merged_into == -1) {
					mark_nodes_as_good(ii);
					break;
				}
			}
		} else {
			for (int ii = 0; ii < clusters.size(); ++ii) {
				if (clusters.get(ii).active) {
					mark_nodes_as_good(ii);
				}
			}
		}

		// if there is no pruning needed, just return the phase 1's result
		if (pruning == PruningMethod.kNoPruning) {
			build_phase1_node_set(phase1_result, result_nodes);
			phase1_result = this.build_phase1_node_set_first_P;
			result_nodes = this.build_phase1_node_set_second_P;
			result_edges = phase1_result;
			return true;
		}
		if (verbosity_level >= 2) {
			System.out.format("------------------------------------------\n");
			System.out.format("Starting pruning\n");
		}
		for (int ii = 0; ii < phase1_result.size(); ++ii) {
			int endPoint0 = edges.get(phase1_result.get(ii))[0];
			int endPoint1 = edges.get(phase1_result.get(ii))[1];
			if (nodeGood.get(endPoint0) && nodeGood.get(endPoint1)) {
				phase2Result.add(phase1_result.get(ii));
			}
		}
		if (pruning == PruningMethod.kSimplePruning) {
			build_phase2_node_set(result_nodes);
			result_nodes = this.build_phase2_node_set_first_P;
			result_edges = phase2Result;
			return true;
		}
		ArrayList<Integer> phase3_result = new ArrayList<Integer>();
		phase3Neighbors = resize(phase3Neighbors, prizes.size(), new ArrayList<Pair<Integer, Double>>());
		for (int ii = 0; ii < phase2Result.size(); ++ii) {
			int cur_edge_index = phase2Result.get(ii);
			int uu = edges.get(cur_edge_index)[0];
			int vv = edges.get(cur_edge_index)[1];
			double cur_cost = costs.get(cur_edge_index);
			phase3Neighbors.get(uu).add(new Pair<Integer, Double>(vv, cur_cost));
			phase3Neighbors.get(vv).add(new Pair<Integer, Double>(uu, cur_cost));
		}

		// GW pruning (this is not strong pruning)
		if (pruning == PruningMethod.kGWPruning) {
			if (verbosity_level >= 2) {
				System.out.format("Starting GW pruning, phase 2 result:\n");
				for (int ii = 0; ii < phase2Result.size(); ++ii) {
					System.out.format("%d ", phase2Result.get(ii));
				}
				System.out.println("\n");
			}
			for (int ii = phase2Result.size() - 1; ii >= 0; --ii) {
				int cur_edge_index = phase2Result.get(ii);
				int uu = edges.get(cur_edge_index)[0];
				int vv = edges.get(cur_edge_index)[1];
				if (nodeDeleted.get(uu) && nodeDeleted.get(vv)) {
					if (verbosity_level >= 2) {
						System.out.format("Not keeping edge %d (%d, %d) because both endpoints already deleted\n",
								cur_edge_index, uu, vv);
					}
					continue;
				}
				if (edgeInfo.get(cur_edge_index).inactiveMergeEvent < 0) {
					mark_clusters_as_necessary(uu);
					mark_clusters_as_necessary(vv);
					phase3_result.add(cur_edge_index);
					if (verbosity_level >= 2) {
						System.out.format("Both endpoint clusters were active, so keeping edge %d (%d, %d)\n",
								cur_edge_index, uu, vv);
					}
				} else {
					InactiveMergeEvent cur_merge_event = inactiveMergeEvents
							.get(edgeInfo.get(cur_edge_index).inactiveMergeEvent);
					int active_side_node = cur_merge_event.active_cluster_node;
					int inactive_side_node = cur_merge_event.inactive_cluster_node;
					int inactive_cluster_index = cur_merge_event.inactive_cluster_index;
					if (clusters.get(inactive_cluster_index).necessary) {
						phase3_result.add(cur_edge_index);
						mark_clusters_as_necessary(inactive_side_node);
						mark_clusters_as_necessary(active_side_node);
						if (verbosity_level >= 2) {
							System.out.format(
									"One endpoint was inactive but is marked necessary (%d), so keeping edge %d (%d, %d)\n",
									inactive_cluster_index, cur_edge_index, uu, vv);
						}
					} else {
						mark_nodes_as_deleted(inactive_side_node, active_side_node);
						if (verbosity_level >= 2) {
							System.out.format(
									"One endpoint was inactive and not marked necessary (%d), so discarding edge %d (%d, %d)\n",
									inactive_cluster_index, cur_edge_index, uu, vv);
						}
					}
				}
			}
			build_phase3_node_set(result_nodes);
			result_nodes = this.build_phase3_node_set_first_P;
			result_edges = phase3_result;
			return true;
			// strong pruning method
		} else if (pruning == PruningMethod.kStrongPruning) {
			if (verbosity_level >= 2) {
				System.out.format("Starting Strong pruning, phase 2 result:\n");
				for (int ii = 0; ii < phase2Result.size(); ++ii) {
					System.out.format("%d ", phase2Result.get(ii));
				}
				System.out.println("\n");
			}
			finalComponentLabel = resize(finalComponentLabel, prizes.size(), -1);
			rootComponentIndex = -1;
			strongPruningParent = resize(strongPruningParent, prizes.size(),
					new Pair<Integer, Double>(-1, -1.0d));
			strongPruningPayoff = resize(strongPruningPayoff, prizes.size(), -1.0);
			for (int ii = 0; ii < phase2Result.size(); ++ii) {
				int cur_node_index = edges.get(phase2Result.get(ii))[0];
				if (finalComponentLabel.get(cur_node_index) == -1) {
					finalComponents.add(new ArrayList<Integer>());
					label_final_component(cur_node_index, finalComponents.size() - 1);
				}
			}
			if (verbosity_level >= 3) {
				System.out.format("number of final_components : %d\n", finalComponents.size());
			}
			for (int ii = 0; ii < (int) finalComponents.size(); ++ii) {
				if (verbosity_level >= 2) {
					System.out.format("Strong pruning on final component %d (size %d):\n", ii,
							finalComponents.get(ii).size());
				}
				if (ii == rootComponentIndex) {
					if (verbosity_level >= 2) {
						System.out.format("Component contains root, pruning starting at %d\n", root);
					}
					strong_pruning_from(root, true);
				} else {
					int best_component_root = find_best_component_root(ii);
					if (verbosity_level >= 3) {
						System.out.format("best_component_root : %d\n", best_component_root);
					}
					if (verbosity_level >= 2) {
						System.out.println("Best start node for current component: " + best_component_root
								+ ", pruning from there\n");
					}
					strong_pruning_from(best_component_root, true);
				}
			}

			for (int ii = 0; ii < phase2Result.size(); ++ii) {
				int cur_edge_index = phase2Result.get(ii);
				int uu = edges.get(cur_edge_index)[0];
				int vv = edges.get(cur_edge_index)[1];
				if (nodeDeleted.get(uu) || nodeDeleted.get(vv)) {
					//////////////////////////////////////////
					if (verbosity_level >= 2) {
						System.out.println("Not keeping edge " + cur_edge_index + " (" + uu + ", " + vv
								+ ") because at least one endpoint already deleted\n");
						// output_function(output_buffer);
					}
					//////////////////////////////////////////
				} else {
					phase3_result.add(cur_edge_index);
				}
			}
			build_phase3_node_set(result_nodes);
			result_nodes = this.build_phase3_node_set_first_P;
			result_edges = phase3_result;
			this.resultEdges = result_edges;
			this.resultNodes = result_nodes;
			return true;
		}
		System.out.println("Error: unknown pruning scheme.\n");
		return false;
	}// run

	public void label_final_component(int start_node_index, int new_component_index) {
		clusterQueue.clear();
		clusterQueue.add(start_node_index);
		finalComponentLabel.set(start_node_index, new_component_index);
		int queue_next = 0;
		while (queue_next < clusterQueue.size()) {
			int cur_node_index = clusterQueue.get(queue_next);
			queue_next += 1;
			finalComponents.get(new_component_index).add(cur_node_index);
			if (cur_node_index == root) {
				rootComponentIndex = new_component_index;
			}
			for (int ii = 0; ii < phase3Neighbors.get(cur_node_index).size(); ++ii) {
				int next_node_index = phase3Neighbors.get(cur_node_index).get(ii).getFirst();
				if (finalComponentLabel.get(next_node_index) == -1) {
					clusterQueue.add(next_node_index);
					finalComponentLabel.set(next_node_index, new_component_index);
				}
			}
		}
	}

	public void strong_pruning_from(int start_node_index, boolean mark_as_deleted) {
		stack.clear();
		stack.add(new Pair<Boolean, Integer>(true, start_node_index));
		strongPruningParent.set(start_node_index, new Pair<Integer, Double>(-1, 0.0));
		if (verbosity_level >= 3) {
			System.out.format("phase3_neighbors size : %d\n", phase3Neighbors.size());
		}
		if (verbosity_level >= 3) {
			System.out.format("stack size is : %d\n", stack.size());
		}
		while (!stack.isEmpty()) {
			int lastElementIndex = stack.size() - 1;
			boolean begin = stack.get(lastElementIndex).getFirst();
			int cur_node_index = stack.get(lastElementIndex).getSecond();
			stack.remove(lastElementIndex);
			if (begin) {
				stack.add(new Pair<Boolean, Integer>(false, cur_node_index));
				for (int ii = 0; ii < phase3Neighbors.get(cur_node_index).size(); ++ii) {
					int next_node_index = phase3Neighbors.get(cur_node_index).get(ii).getFirst();
					double next_cost = phase3Neighbors.get(cur_node_index).get(ii).getSecond();
					if (next_node_index == strongPruningParent.get(cur_node_index).getFirst()) {
						continue;
					}
					strongPruningParent.set(next_node_index, new Pair<Integer, Double>(cur_node_index, next_cost));
					stack.add(new Pair<Boolean, Integer>(true, next_node_index));
				}
			} else {
				strongPruningPayoff.set(cur_node_index, prizes.get(cur_node_index));
				for (int ii = 0; ii < phase3Neighbors.get(cur_node_index).size(); ++ii) {
					int next_node_index = phase3Neighbors.get(cur_node_index).get(ii).getFirst();
					double next_cost = phase3Neighbors.get(cur_node_index).get(ii).getSecond();
					if (next_node_index == strongPruningParent.get(cur_node_index).getFirst().intValue()) {
						continue;
					}
					double next_payoff = strongPruningPayoff.get(next_node_index) - next_cost;
					if (next_payoff <= 0.0) {
						if (mark_as_deleted) {
							if (verbosity_level >= 2) {
								System.out.format(
										"Subtree starting at %d has a nonpositive contribution of %6f, pruning (good side: %d)\n",
										next_node_index, next_payoff, cur_node_index);
							}
							mark_nodes_as_deleted(next_node_index, cur_node_index);
						}
					} else {
						double currentValue = strongPruningPayoff.get(cur_node_index);
						strongPruningPayoff.set(cur_node_index, currentValue + next_payoff);
					}
				}
			}
		}
	}

	public int find_best_component_root(int component_index) {
		int cur_best_root_index = finalComponents.get(component_index).get(0);
		strong_pruning_from(cur_best_root_index, false);
		double cur_best_value = strongPruningPayoff.get(cur_best_root_index);
		stack2.clear();
		for (int ii = 0; ii < phase3Neighbors.get(cur_best_root_index).size(); ++ii) {
			stack2.add(phase3Neighbors.get(cur_best_root_index).get(ii).getFirst());
		}
		while (!stack2.isEmpty()) {
			int cur_node_index = stack2.get(stack2.size() - 1);
			stack2.remove(stack2.size() - 1);
			int cur_parent_index = strongPruningParent.get(cur_node_index).getFirst();
			double parent_edge_cost = strongPruningParent.get(cur_node_index).getSecond();
			double parent_val_without_cur_node = strongPruningPayoff.get(cur_parent_index);
			double cur_node_net_payoff = strongPruningPayoff.get(cur_node_index) - parent_edge_cost;
			if (cur_node_net_payoff > 0.0) {
				parent_val_without_cur_node -= cur_node_net_payoff;
			}
			if (parent_val_without_cur_node > parent_edge_cost) {
				double currentPayOff = strongPruningPayoff.get(cur_node_index);
				strongPruningPayoff.set(cur_node_index,
						currentPayOff + (parent_val_without_cur_node - parent_edge_cost));
			}
			if (strongPruningPayoff.get(cur_node_index) > cur_best_value) {
				cur_best_root_index = cur_node_index;
				cur_best_value = strongPruningPayoff.get(cur_node_index);
			}
			for (int ii = 0; ii < phase3Neighbors.get(cur_node_index).size(); ++ii) {
				int next_node_index = phase3Neighbors.get(cur_node_index).get(ii).getFirst();
				if (next_node_index != cur_parent_index) {
					stack2.add(next_node_index);
				}
			}
		}
		return cur_best_root_index;
	}

	public void build_phase3_node_set(ArrayList<Integer> node_set) {
		node_set.clear();
		for (int ii = 0; ii < (int) prizes.size(); ++ii) {
			if (!nodeDeleted.get(ii) && nodeGood.get(ii)) {
				node_set.add(ii);
			}
		}
		this.build_phase3_node_set_first_P = node_set;
	}

	public void build_phase2_node_set(ArrayList<Integer> node_set) {
		node_set.clear();
		for (int ii = 0; ii < (int) prizes.size(); ++ii) {
			if (nodeGood.get(ii)) {
				node_set.add(ii);
			}
		}
		this.build_phase2_node_set_first_P = node_set;
	}

	public void build_phase1_node_set(ArrayList<Integer> edge_set, ArrayList<Integer> node_set) {
		ArrayList<Boolean> included = new ArrayList<Boolean>();
		for (int i = 0; i < prizes.size(); i++) {
			included.add(false);
		}
		node_set.clear();
		for (int ii = 0; ii < edge_set.size(); ++ii) {
			int uu = edges.get(edge_set.get(ii))[0];
			int vv = edges.get(edge_set.get(ii))[1];

			if (!included.get(uu)) {
				included.set(uu, true);
				node_set.add(uu);
			}
			if (!included.get(vv)) {
				included.set(vv, true);
				node_set.add(vv);
			}
		}
		for (int ii = 0; ii < prizes.size(); ++ii) {
			if (nodeGood.get(ii) && (!included.get(ii))) {
				node_set.add(ii);
			}
		}
		this.build_phase1_node_set_first_P = edge_set;
		this.build_phase1_node_set_second_P = node_set;
	}

	public void get_statistics(FastPCST.Statistics s) {
		s = stats;
	}

	/** PruningMethod */
	public enum PruningMethod {
		kNoPruning(0), kSimplePruning(1), kGWPruning(2), kStrongPruning(3), kUnknownPruning(4);

		private int intValue;
		private static HashMap<Integer, PruningMethod> mappings;

		private static HashMap<Integer, PruningMethod> getMappings() {
			if (mappings == null) {
				synchronized (PruningMethod.class) {
					if (mappings == null) {
						mappings = new HashMap<Integer, PruningMethod>();
					}
				}
			}
			return mappings;
		}

		private PruningMethod(int value) {
			intValue = value;
			getMappings().put(value, this);
		}

		public int getValue() {
			return intValue;
		}

		public static PruningMethod forValue(int value) {
			return getMappings().get(value);
		}
	}
	
	
	private ArrayList<EdgeInfo> resize(ArrayList<EdgeInfo> arr, int size, EdgeInfo val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<EdgeInfo>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<EdgeInfo>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(new EdgeInfo());
            }
            return arr;
        }
        return arr;
    }

	private ArrayList<Integer> resize(ArrayList<Integer> arr, int size, Integer val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<Integer>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<Integer>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(val);
            }
            return arr;
        }
        return arr;
    }

    private ArrayList<Boolean> resize(ArrayList<Boolean> arr, int size, Boolean val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<Boolean>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<Boolean>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(val);
            }
            return arr;
        }
        return arr;
    }


    private ArrayList<Pair<Integer, Double>> resize(
            ArrayList<Pair<Integer, Double>> arr, int size, Pair<Integer, Double> pair) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<Pair<Integer, Double>>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<Pair<Integer, Double>>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(pair);
            }
            return arr;
        }
        return arr;
    }

    private ArrayList<ArrayList<Pair<Integer, Double>>> resize(
            ArrayList<ArrayList<Pair<Integer, Double>>> arr,
            int size, ArrayList<Pair<Integer, Double>> pair) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<ArrayList<Pair<Integer, Double>>>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<ArrayList<Pair<Integer, Double>>>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(new ArrayList<Pair<Integer, Double>>());
            }
            return arr;
        }
        return arr;
    }

    private ArrayList<Double> resize(ArrayList<Double> arr, int size, Double val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<Double>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<Double>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(val);
            }
            return arr;
        }
        return arr;
    }

    private ArrayList<EdgePart> resize(ArrayList<EdgePart> arr, int size, EdgePart val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<EdgePart>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<EdgePart>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(new EdgePart());
            }
            return arr;
        }
        return arr;
    }

	public class Statistics {
		public long total_num_edge_events;
		public long num_deleted_edge_events;
		public long num_merged_edge_events;
		public long total_num_merge_events;
		public long num_active_active_merge_events;
		public long num_active_inactive_merge_events;
		public long total_num_edge_growth_events;
		public long num_active_active_edge_growth_events;
		public long num_active_inactive_edge_growth_events;
		public long num_cluster_events;
	}

	public class InactiveMergeEvent {
		public int active_cluster_index;
		public int inactive_cluster_index;
		public int active_cluster_node;
		public int inactive_cluster_node;
	}

	public class Cluster {
		public PairingHeap edge_parts;
		public boolean active;
		public double active_start_time;
		public double active_end_time;
		public int merged_into;
		public double prize_sum;
		public double subcluster_moat_sum;
		public double moat;
		public boolean contains_root;
		public int skip_up;
		public double skip_up_sum;
		public int merged_along;
		public int child_cluster_1;
		public int child_cluster_2;
		public boolean necessary;

		public Cluster(ArrayList<PairingHeap.Node> heap_buffer) {
			this.edge_parts = new PairingHeap(heap_buffer);
		}
	}
}
