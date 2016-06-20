package edu.albany.cs.base;

import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.Edge;
import edu.albany.cs.base.Utils;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * In this class, using APMDIOFormat class/object to get the data from like-stp
 * file. How to use it ? Specify your input file
 *
 * @author Baojian Zhou
 */
public class APDMInputFormat {

	// header
	private final static String inputHeader = Utils.commentLine + "\n"
			+ "#APDM Input Graph, this input graph includes 3 sections:\n" + "#section1 : general information\n"
			+ "#section2 : nodes\n" + "#section3 : edges\n" + "#section4 : trueSubGraph (Optional)\n" + "#\n"
			+ "#if nodes haven't information set weight to null\n" + Utils.commentLine + "\n";

	public InputData data;

	/**
	 * read APDM input file
	 *
	 * @param APDMFile
	 *            APDM input file
	 */
	public APDMInputFormat(File APDMFile) {

		/** read inputFile */
		readAPDMFile(APDMFile);

		/** initialize true subgraph nodes */
		if (data.trueSubGraphEdges != null) {
			int[] trueSubGraphNodes = null;
			for (int[] edge : data.trueSubGraphEdges.keySet()) {
				if (!ArrayUtils.contains(trueSubGraphNodes, edge[0])) {
					trueSubGraphNodes = ArrayUtils.add(trueSubGraphNodes, edge[0]);
				}
				if (!ArrayUtils.contains(trueSubGraphNodes, edge[1])) {
					trueSubGraphNodes = ArrayUtils.add(trueSubGraphNodes, edge[1]);
				}
			}
			data.trueSubGraphNodes = trueSubGraphNodes;
		} else {
			data.trueSubGraphNodes = null;
		}

		/** initialize adjacency information */
		ArrayList<ArrayList<Integer>> graphAdjList = new ArrayList<>();
		ArrayList<ArrayList<Double>> graphWeightedAdjList = new ArrayList<>();
		int[][] graphAdj = new int[data.numNodes][];
		double[][] graphWeightedAdj = new double[data.numNodes][];
		for (int i = 0; i < data.numNodes; i++) {
			graphAdjList.add(new ArrayList<Integer>());
			graphWeightedAdjList.add(new ArrayList<Double>());
		}

		for (int[] edge : data.edges.keySet()) {
			if (!graphAdjList.get(edge[0]).contains(edge[1])) {
				graphAdjList.get(edge[0]).add(edge[1]);
				graphWeightedAdjList.get(edge[0]).add(data.edges.get(edge));
				graphAdj[edge[0]] = ArrayUtils.add(graphAdj[edge[0]], edge[1]);
				graphWeightedAdj[edge[0]] = ArrayUtils.add(graphWeightedAdj[edge[0]], data.edges.get(edge));
			}
			if (!graphAdjList.get(edge[1]).contains(edge[0])) {
				graphAdjList.get(edge[1]).add(edge[0]);
				graphWeightedAdjList.get(edge[1]).add(data.edges.get(edge));
				graphAdj[edge[1]] = ArrayUtils.add(graphAdj[edge[1]], edge[0]);
				graphWeightedAdj[edge[1]] = ArrayUtils.add(graphWeightedAdj[edge[1]], data.edges.get(edge));
			}
		}
		data.graphAdj = graphAdj;
		data.graphAdjList = graphAdjList;
		data.graphWeightedAdj = graphWeightedAdj;
		data.graphWeightedAdjList = graphWeightedAdjList;

		/** other parameters. */
		data.V = new int[data.numNodes];
		for (int i = 0; i < data.V.length; i++) {
			data.V[i] = i;
		}
		data.identityb = new double[data.numNodes];
		Arrays.fill(data.identityb, 1.0D);
		data.smallb = new double[data.numNodes];
		Arrays.fill(data.smallb, 0.01D);
		data.cc = new ConnectedComponents(data.graphAdjList);

	}

	public APDMInputFormat(String apdmFileName) {
		this(new File(apdmFileName));/** read inputFile */
	}

	/**
	 * According to the APDM input format, we generate APDM input file
	 *
	 * @param usedAlgorithm
	 *            this kind of input file is used for usedAlgorithm
	 * @param dataSource
	 *            dataSource [WaterPollutionDataset CivilUnrestDataset
	 *            GridDataset SNAPDataset TransportationDataset ]
	 * @param edges
	 *            the edges of the graph
	 * @param PValue
	 *            the nodes of the graph
	 * @param fileName
	 *            the input file name
	 */
	public static void generateAPDMFile(String usedAlgorithm, String dataSource, ArrayList<Edge> edges, double[] PValue,
			double[] counts, double[] averValue, HashMap<int[], Double> trueSubGraphEdges, String fileName) {
		DecimalFormat decimalFormat = new DecimalFormat("0.000000");
		try {
			FileWriter fw;
			fw = new FileWriter(fileName, false);
			fw.write(APDMInputFormat.inputHeader);
			// general information
			fw.write("SECTION1 (General Information)\n");
			if (PValue == null) {
				fw.write("numNodes = " + 0 + "\n");
			} else {
				fw.write("numNodes = " + PValue.length + "\n");
			}
			if (edges == null) {
				fw.write("numEdges = " + 0 + "\n");
			} else {
				fw.write("numEdges = " + edges.size() + "\n");
			}
			fw.write("usedAlgorithm = " + usedAlgorithm + "\n");
			fw.write("dataSource = " + dataSource + "\n");
			fw.write("END\n" + Utils.commentLine + "\n");

			// nodes information
			fw.write("SECTION2 (Nodes Information)\n");
			fw.write("NodeID PValue Counts\n");
			if (PValue == null) {
				fw.write("null\n");
			} else {
				if (counts == null && averValue == null) {
					for (int i = 0; i < PValue.length; i++) {
						fw.write(i + " " + decimalFormat.format(PValue[i]) + "\n");
					}
				} else if (counts != null && averValue == null) {
					for (int i = 0; i < PValue.length; i++) {
						fw.write(i + " " + decimalFormat.format(PValue[i]) + " " + counts[i] + "\n");
					}
				} else {
					for (int i = 0; i < PValue.length; i++) {
						fw.write(i + " " + decimalFormat.format(PValue[i]) + " " + counts[i] + " " + averValue[i]
								+ "\n");
					}
				}

			}
			fw.write("END\n" + Utils.commentLine + "\n");

			// edges information
			fw.write("SECTION3 (Edges Information)\n");
			fw.write("EndPoint0 EndPoint1 Weight\n");
			if (edges != null) {
				for (Edge e : edges) {
					fw.write(e.i + " " + e.j + " " + decimalFormat.format(e.cost) + "\n");
				}
			}
			fw.write("END\n" + Utils.commentLine + "\n");
			// edges information
			fw.write("SECTION4 (TrueSubGraph Information)\n");
			fw.write("EndPoint0 EndPoint1 Weight\n");
			if (trueSubGraphEdges != null) {
				for (int[] e : trueSubGraphEdges.keySet()) {
					fw.write(e[0] + " " + e[1] + " " + trueSubGraphEdges.get(e) + "\n");
				}
			}
			fw.write("END\n" + Utils.commentLine + "\n");
			fw.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
	}

	/**
	 * According to the APDM input format, we generate APDM input file
	 *
	 * @param usedAlgorithm
	 *            this kind of input file is used for usedAlgorithm
	 * @param dataSource
	 *            dataSource [WaterPollutionDataset CivilUnrestDataset
	 *            GridDataset SNAPDataset TransportationDataset ]
	 * @param edges
	 *            the edges of the graph
	 * @param PValue
	 *            the nodes of the graph
	 * @param fileName
	 *            the input file name
	 */
	public static void generateAPDMFile(String usedAlgorithm, String dataSource, ArrayList<Edge> edges, double[] PValue,
			double[] counts, double[] averValue, ArrayList<Edge> trueSubGraphEdges, String fileName)
					throws IOException {
		DecimalFormat decimalFormat = new DecimalFormat("0.000000");

		FileWriter fw = new FileWriter(fileName, false);
		fw.write(APDMInputFormat.inputHeader);
		// general information
		fw.write("SECTION1 (General Information)\n");
		if (PValue == null) {
			fw.write("numNodes = " + 0 + "\n");
		} else {
			fw.write("numNodes = " + PValue.length + "\n");
		}
		if (edges == null) {
			fw.write("numEdges = " + 0 + "\n");
		} else {
			fw.write("numEdges = " + edges.size() + "\n");
		}
		fw.write("usedAlgorithm = " + usedAlgorithm + "\n");
		fw.write("dataSource = " + dataSource + "\n");
		fw.write("END\n" + Utils.commentLine + "\n");

		// nodes information
		fw.write("SECTION2 (Nodes Information)\n");
		fw.write("NodeID PValue Counts\n");
		if (PValue == null) {
			fw.write("null\n");
		} else {
			if (counts == null && averValue == null) {
				for (int i = 0; i < PValue.length; i++) {
					fw.write(i + " " + decimalFormat.format(PValue[i]) + "\n");
				}
			} else {
				for (int i = 0; i < PValue.length; i++) {
					fw.write(i + " " + decimalFormat.format(PValue[i]) + " " + counts[i] + " " + averValue[i] + "\n");
				}
			}

		}
		fw.write("END\n" + Utils.commentLine + "\n");

		// edges information
		fw.write("SECTION3 (Edges Information)\n");
		fw.write("EndPoint0 EndPoint1 Weight\n");
		if (edges != null) {
			for (Edge e : edges) {
				fw.write(e.i + " " + e.j + " " + decimalFormat.format(e.cost) + "\n");
			}
		}
		fw.write("END\n" + Utils.commentLine + "\n");
		// edges information
		fw.write("SECTION4 (TrueSubGraph Information)\n");
		fw.write("EndPoint0 EndPoint1 Weight\n");
		if (trueSubGraphEdges != null) {
			for (Edge e : trueSubGraphEdges) {
				fw.write(e.i + " " + e.j + " " + e.cost + "\n");
			}
		}
		fw.write("END\n" + Utils.commentLine + "\n");
		fw.close();
	}

	/**
	 * @param APDMFile
	 *            read APDM from file
	 */
	private boolean readAPDMFile(File APDMFile) {
		BufferedReader br = null;
		data = new InputData();
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(APDMFile));

			while ((sCurrentLine = br.readLine()) != null) {
				if (sCurrentLine.startsWith("#")) {
					continue;
				}
				// general Information
				if (sCurrentLine.startsWith("SECTION1")) {
					while (!(sCurrentLine = br.readLine()).equals("END")) {
						if (sCurrentLine.startsWith("numNodes")) {
							String[] str = sCurrentLine.split(" ");
							data.numNodes = Integer.parseInt(str[2]);
						}
						if (sCurrentLine.startsWith("numEdges")) {
							String[] str = sCurrentLine.split(" ");
							int numEdges = Integer.parseInt(str[2]);
							data.numEdges = numEdges;
						}
						if (sCurrentLine.startsWith("dataSource")) {
							String[] str = sCurrentLine.trim().split(" ");
							String dataSource = str[2];
							this.data.dataSource = dataSource;
						}
						if (sCurrentLine.startsWith("usedAlgorithm")) {
							String[] str = sCurrentLine.split(" ");
							String usedAlgorithm = str[2];
							this.data.usedAlgorithm = usedAlgorithm;
						}
					}
				}

				// nodes information
				if (sCurrentLine.startsWith("SECTION2")) {
					data.nodes = new HashMap<Integer, Double>();
					data.PValue = new double[data.numNodes];
					data.counts = new double[data.numNodes];
					data.averCounts = new double[data.numNodes];
					data.speed = new double[data.numNodes];
					data.meanSpeed = new double[data.numNodes];
					data.std = new double[data.numNodes];
					data.mean = new double[data.numNodes];
					int count = 0;
					while (!(sCurrentLine = br.readLine()).equals("END")) {
						if (sCurrentLine.startsWith("NodeID")) {
							continue;
						}
						String[] str = sCurrentLine.split(" ");
						data.nodes.put(Integer.parseInt(str[0]), Double.parseDouble(str[1]));
						data.PValue[count] = Double.parseDouble(str[1]);
						if (data.dataSource.equals("CivilUnrest") || data.dataSource.equals("citHepPh")) {
							data.counts[count] = Double.parseDouble(str[2]);
							if (data.dataSource.equals("citHepPh")) {
								data.averCounts[count] = Double.parseDouble(str[3]);
							}
						} else if (data.dataSource.equals("Transportation")) {
							data.speed[count] = Double.parseDouble(str[2]);
							data.meanSpeed[count] = Double.parseDouble(str[3]);
						} else if (data.dataSource.equals("Haze") || data.dataSource.equals("CrimeOfChicago")) {
							data.counts[count] = Double.parseDouble(str[2]);
							data.averCounts[count] = Double.parseDouble(str[3]);
						} else if (data.dataSource.equals("WaterDataSet") || data.dataSource.equals("WaterData")) {
							data.counts[count] = Double.parseDouble(str[2]);
						} else if (data.dataSource.equals("DiseaseOutBreakSimu")) {
							data.counts[count] = Double.parseDouble(str[1]);
							data.mean[count] = Double.parseDouble(str[2]);
							data.std[count] = Double.parseDouble(str[3]);
						} else if (data.dataSource.equals("NewYorkCityTaxi")) {
							data.counts[count] = Double.parseDouble(str[1]);
						}
						count++;
					}
				}
				// edges information
				if (sCurrentLine.startsWith("SECTION3")) {
					data.edges = new HashMap<int[], Double>();
					data.newEdges = new ArrayList<Edge>();
					data.intEdges = new ArrayList<Integer[]>();
					data.edgeCosts = new ArrayList<Double>();
					int count = 0;
					while (!(sCurrentLine = br.readLine()).equals("END")) {
						if (sCurrentLine.startsWith("EndPoint0")) {
							continue;
						}
						String[] str = sCurrentLine.split(" ");
						int[] edge = new int[] { Integer.parseInt(str[0]), Integer.parseInt(str[1]) };
						data.edges.put(edge, Double.parseDouble(str[2]));
						data.newEdges.add(new Edge(edge[0], edge[1], count++, Double.parseDouble(str[2])));
						data.intEdges.add(new Integer[] { edge[0], edge[1] });
						data.edgeCosts.add(Double.parseDouble(str[2]));
					}
					ArrayList<Double> identityEdgeCosts = new ArrayList<Double>();
					for (int i = 0; i < data.edgeCosts.size(); i++) {
						identityEdgeCosts.add(1.0D);
					}
					data.identityEdgeCosts = identityEdgeCosts;
				}
				// trueSubGraphEdges information
				if (sCurrentLine.startsWith("SECTION4")) {
					data.trueSubGraphEdges = new HashMap<int[], Double>();
					while (!(sCurrentLine = br.readLine()).equals("END")) {
						if (sCurrentLine.startsWith("EndPoint0") || sCurrentLine.startsWith("null")) {
							continue;
						}
						String[] str = sCurrentLine.split(" ");
						int[] edge = new int[] { Integer.parseInt(str[0]), Integer.parseInt(str[1]) };
						data.trueSubGraphEdges.put(edge, Double.parseDouble(str[2]));
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return true;
	}

	// test
	public static void main(String args[]) throws IOException {
		processCrimeOfChicago();
	}

	public static void processCrimeOfChicago() throws IOException {

		for (File file : new File("./data/CrimeOfChicago/graph").listFiles()) {
			APDMInputFormat apdm = new APDMInputFormat(file);
			double[] c = apdm.data.counts;
			System.out.println(StatUtils.sum(c));
			double mu = StatUtils.mean(c);
			double sigma = Math.sqrt(StatUtils.variance(c));
			System.out.println(mu + " " + sigma);
			double[] PValue = new double[c.length];
			NormalDistribution normal = new NormalDistribution(mu, sigma);
			for (int i = 0; i < c.length; i++) {
				PValue[i] = 1 - normal.cumulativeProbability(c[i]);
				if (PValue[i] <= 0.0D) {
					PValue[i] = 0.0001D;
				}
			}
			APDMInputFormat.generateAPDMFile("NULL", apdm.data.dataSource, apdm.data.newEdges, PValue, apdm.data.counts,
					apdm.data.averCounts, new HashMap<int[], Double>(), "./tmp/" + file.getName());
		}
	}

	public static void processCitHepPh(String fileName, String year) throws IOException {
		int[] nodes = Utils.getIntFromFile("./data/citHepPh/maximum_CC_nodes.txt");
		ArrayList<Integer> nodesID = new ArrayList<Integer>();
		for (int i : nodes) {
			nodesID.add(i);
		}
		HashMap<Integer, Double> pvalue = new HashMap<Integer, Double>();
		HashMap<Integer, Integer> citCnts = new HashMap<Integer, Integer>();
		pvalue = Utils.getPValueFromFile("./data/citHepPh/" + fileName);
		citCnts = Utils.getCitCntsFromFile("./data/citHepPh/" + fileName);
		double[] PValue = new double[nodesID.size()];
		double[] counts = new double[nodesID.size()];
		double[] average = new double[nodesID.size()];
		for (int i = 0; i < average.length; i++) {
			average[i] = 2.391845;
		}
		for (int i = 0; i < PValue.length; i++) {
			PValue[i] = pvalue.get(nodesID.get(i));
			counts[i] = citCnts.get(nodesID.get(i));
		}
		ArrayList<Integer[]> edges = Utils.getGraphEdgeFromFile("data/citHepPh/edges.txt", " ");
		ArrayList<Edge> Edges = new ArrayList<Edge>();
		int count = 0;
		for (Integer[] edge : edges) {
			Edges.add(new Edge(nodesID.indexOf(edge[0]), nodesID.indexOf(edge[1]), count++, 1.0D));
		}
		// APDMInputFormat.generateAPDMFile("NULL", "citHepPh", Edges, PValue,
		// "./realDataSet/citHepPh/APDM-Graph-CitHepPh-2002-11895.txt") ;
		APDMInputFormat.generateAPDMFile("NULL", "citHepPh", Edges, PValue, counts, average,
				new HashMap<int[], Double>(), "./data/citHepPh/APDM-Graph-CitHepPh-" + year + "-11895.txt");
	}

}
