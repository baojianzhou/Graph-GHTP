package edu.albany.cs.base;

import org.apache.commons.lang3.ArrayUtils;

import edu.albany.cs.scoreFuncs.FuncType;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * this Utils class is for testing and printing some results
 *
 * @author Baojian
 */
public class Utils {
	public static int count = 0;

	public static void testLine() {
		System.out.println("test line " + count + " ...");
		count++;
	}

	public final static String Line = "\n=======================================================\n";
	public final static String InputLine = "\n========================= Input  ==============================\n";
	public final static String OutputLine = "\n========================= Output ==============================\n";
	public final static String commentLine = "#################################################################";

	public static long getPID() {
		String processName = java.lang.management.ManagementFactory.getRuntimeMXBean().getName();
		return Long.parseLong(processName.split("@")[0]);
	}

	public static double median(double[] array) {
		Arrays.sort(array);
		double median;
		if (array.length % 2 == 0)
			median = (array[array.length / 2] + array[array.length / 2 - 1]) / 2;
		else
			median = array[array.length / 2];
		return median;
	}

	public static String calPreAndRec(int[] resultNodes, int[] trueNodes) {
		int[] intersect = Utils.intersect(resultNodes, trueNodes);
		if (intersect == null || intersect.length == 0) {
			System.out.println("precision : " + intersect.length * 1.0D / resultNodes.length * 1.0D);
			System.out.println("recall : " + intersect.length * 1.0D / trueNodes.length * 1.0D);
		} else {
			System.out.println("precision : " + intersect.length * 1.0D / resultNodes.length * 1.0D);
			System.out.println("recall : " + intersect.length * 1.0D / trueNodes.length * 1.0D);
		}
		double precision = (intersect.length * 1.0D / resultNodes.length * 1.0D);
		double recall = (intersect.length * 1.0D / trueNodes.length * 1.0D);
		DecimalFormat df = new DecimalFormat("#.###");

		String str = "precision : " + df.format(precision) + " recall : " + df.format(recall) + "  ====  ";

		return str;
	}

	public static double[] calPreAndRecDouble(int[] resultNodes, int[] trueNodes) {
		int[] intersect = Utils.intersect(resultNodes, trueNodes);
		double[] pre_rec = new double[2];
		if (intersect == null || intersect.length == 0) {
			System.out.println("precision : " + intersect.length * 1.0D / resultNodes.length * 1.0D);
			System.out.println("recall : " + intersect.length * 1.0D / trueNodes.length * 1.0D);
		} else {
			System.out.println("precision : " + intersect.length * 1.0D / resultNodes.length * 1.0D);
			System.out.println("recall : " + intersect.length * 1.0D / trueNodes.length * 1.0D);
		}
		pre_rec[0] = (intersect.length * 1.0D / resultNodes.length * 1.0D);
		pre_rec[1] = (intersect.length * 1.0D / trueNodes.length * 1.0D);
		return pre_rec;
	}

	public static void writeResultData(String outputFileName, double[] c, int[] trueSubGraphNodes,
			int[] resultNodes_supportX, double pre, double rec, double funcValue) throws IOException {

		FileWriter fileWriter = new FileWriter(outputFileName, false);

		// The first line is all abnormal nodes (including noise nodes)
		String str = "";
		for (int i = 0; i < c.length; i++) {
			if (c[i] != 0.0D) {
				str += i + " ";
			}
		}
		str += "\n";
		fileWriter.write(str);

		// The second line is true abnormal nodes
		str = "";
		for (int i : trueSubGraphNodes) {
			str += i + " ";
		}
		str += "\n";
		fileWriter.write(str);
		// The third line is nodes of result
		str = "";
		for (int k : resultNodes_supportX) {
			str += k + " ";
		}
		str += "\n";
		fileWriter.write(str);
		// The fourth line is edges in forest
		for (ArrayList<Integer[]> tree : new ArrayList<ArrayList<Integer[]>>()) {
			str = "";
			for (Integer[] edge : tree) {
				str += edge[0] + " " + edge[1] + " ;";
			}
			fileWriter.write(str);
			fileWriter.write("#");
		}
		fileWriter.write("corresponding recall : " + rec + " ; precision : " + pre + "\n");
		fileWriter.write("function value : " + funcValue + "\n");
		fileWriter.close();
	}

	public static int[] getIntArrayFromIntegerList(ArrayList<Integer> list) {

		if (list == null) {
			return null;
		}
		if (list.isEmpty()) {
			return new int[0];
		}
		int[] resultArray = new int[list.size()];
		for (int i = 0; i < list.size(); i++) {
			resultArray[i] = list.get(i);
		}
		return resultArray;
	}

	public static String getOsName() {
		String OS = null;
		if (OS == null) {
			OS = System.getProperty("os.name");
		}
		return OS;
	}

	public static void writeFile(String fileName, String content) {
		try {
			FileWriter fw = new FileWriter(fileName, true);
			fw.write(content + "\n");
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static boolean isWindows() {
		return getOsName().startsWith("Windows");
	}

	public static boolean isLinux() {
		return getOsName().startsWith("Linux");
	}

	/**
	 * print Integer array
	 *
	 * @param A
	 *            integer array
	 */
	public static void print(Integer[] A) {
		int count = 20;
		System.out.println();
		if (A != null) {
			for (int i = 0; i < A.length; i++) {
				System.out.print(A[i] + " ");
				if (i % count == 19) {
					System.out.println();
				}
			}
		} else {
			System.out.println("\nthe set is empty");
		}
	}
	
	/**
	 * print Integer array
	 *
	 * @param A
	 *            integer array
	 */
	public static void print(ArrayList<Integer> A) {
		int count = 20;
		System.out.println();
		if (A != null) {
			for (int i = 0; i < A.size(); i++) {
				System.out.print(A.get(i) + " ");
				if (i % count == 19) {
					System.out.println();
				}
			}
		} else {
			System.out.println("\nthe set is empty");
		}
	}

	/**
	 * print int array (this array should be less than 1000)
	 *
	 * @param A
	 *            integer array
	 */
	public static void print(int[] A) {
		int count = 20;
		System.out.println();
		if (A != null) {
			if (A.length <= 4000) {
				for (int i = 0; i < A.length; i++) {
					System.out.print(A[i] + " ");
					if (i % count == 0 && i != 0) {
						System.out.println();
					}
				}
			} else {
				System.out.println("the length of array is larger than 1000, don' show it.");
			}
		} else {
			System.out.println("\nthe set is empty");
		}
		System.out.println();
	}

	/**
	 * print int array, one line has count number of int value
	 *
	 * @param A
	 *            integer array
	 * @param count
	 *            number of count
	 */
	public static void print(int[] A, int count) {
		if (A != null) {
			for (int i = 0; i < A.length; i++) {
				System.out.print(A[i] + " ");
				if (i % count == 0) {
					System.out.println();
				}
			}
		} else {
			System.out.println("\nthe set is empty");
		}
	}

	/**
	 * print double array
	 *
	 * @param A
	 *            integer array
	 */
	public static void print(double[] A) {
		int count = 20;
		System.out.println();
		for (int i = 0; i < A.length; i++) {
			System.out.print(A[i] + " ");
			if (i % count == 19) {
				System.out.println();
			}
		}
	}

	/**
	 * print maxtrix
	 *
	 * @param maxtrix
	 *            matrix
	 */
	public static void printMatrix(Integer[][] maxtrix) {
		if (maxtrix != null) {
			for (int i = 0; i < maxtrix.length; i++) {
				for (int j = 0; j < maxtrix.length; j++) {
					System.out.print(maxtrix[i][j] + " ");
				}
				System.out.println();
			}
		} else {
			System.out.println("this is empty maxtrix");
		}
	}

	/**
	 * print maxtrix
	 */
	public static void printMatrix(double[][] maxtrix) {
		if (maxtrix != null) {
			for (int i = 0; i < maxtrix.length; i++) {
				for (int j = 0; j < maxtrix.length; j++) {
					System.out.print(maxtrix[i][j] + " ");
				}
				System.out.println();
			}
		} else {
			System.out.println("this is empty maxtrix");
		}
	}

	/**
	 * return the intersection of a and b
	 *
	 * @return a intersect b
	 */
	public static Integer[] intersect(Integer[] a, Integer[] b) {
		Integer[] result = new Integer[Math.min(a.length, b.length)];
		HashSet<Integer> hs = new HashSet<Integer>();
		for (int i = 0; i < b.length; i++) {
			hs.add(b[i]);
		}
		int count = 0;
		for (int i = 0; i < a.length; i++) {
			if (hs.contains(a[i])) {
				result[count++] = a[i];
			}
		}
		return ArrayUtils.subarray(result, 0, count);
	}

	/**
	 * return the intersection of a and b
	 *
	 * @return a intersect b
	 */
	public static int[] intersect(int[] a, int[] b) {
		if (a == null || b == null) {
			return null;
		}
		int[] result = new int[Math.min(a.length, b.length)];
		HashSet<Integer> hs = new HashSet<Integer>();
		for (int i = 0; i < b.length; i++) {
			hs.add(b[i]);
		}
		int count = 0;
		for (int i = 0; i < a.length; i++) {
			if (hs.contains(a[i])) {
				result[count++] = a[i];
			}
		}
		return ArrayUtils.subarray(result, 0, count);
	}

	/**
	 * return the intersection of a and b
	 *
	 * @return a intersect b
	 */
	public static Integer[] intersect(Integer[] a, int[] b) {
		if (a == null || b == null) {
			return null;
		} else {
			Integer[] result = new Integer[Math.min(a.length, b.length)];
			HashSet<Integer> hs = new HashSet<Integer>();
			for (int i = 0; i < b.length; i++) {
				hs.add(b[i]);
			}
			int count = 0;
			for (int i = 0; i < a.length; i++) {
				if (hs.contains(a[i])) {
					result[count++] = a[i];
				}
			}
			return ArrayUtils.subarray(result, 0, count);
		}
	}

	/**
	 * return the difference of a and b, Integer array
	 *
	 * @return a \ b
	 */
	public static Integer[] setADifferentB(Integer[] a, Integer[] b) {
		Integer[] result = null;
		for (int i : a) {
			if (!ArrayUtils.contains(b, i)) {
				result = ArrayUtils.add(result, i);
			}
		}
		return result;
	}

	/**
	 * return the difference of a and b, int array
	 */
	public static int[] setADifferentB(int[] a, int[] b) {
		int[] result = null;
		for (int i : a) {
			if (!ArrayUtils.contains(b, i)) {
				result = ArrayUtils.add(result, i);
			}
		}
		return result;
	}

	/**
	 * stop for testing
	 */
	public static void stop() {
		try {
			System.out.println("Press any key to continue ...");
			System.in.read();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void delay() {
		try {
			Thread.sleep(3000);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	/**
	 * print the ArrayList, adjList stands for a graph
	 */
	public static void printArrayList(ArrayList<ArrayList<Integer>> adjList) {
		if (adjList != null) {
			for (int i = 0; i < adjList.size(); i++) {
				for (int k : adjList.get(i)) {
					System.out.print(k + ", ");
				}
				System.out.println();
			}
		} else {
			System.out.println("this is an empty adjList");
		}
	}

	/**
	 * for the shortest path, we print the path
	 */
	public static void print(HashSet<Integer> path) {
		Iterator<Integer> iter = path.iterator();
		int count = 0;
		while (iter.hasNext()) {
			System.out.println(iter.next());
			if ((++count) % 50 == 0)
				System.out.println();
		}
	}

	/**
	 * get all the nodes in the fileName, it will return Integer array
	 *
	 * @return all the nodes in the fileName
	 */
	public static Integer[] getIntegerNodesFromFile(String fileName) {
		BufferedReader br = null;
		Integer[] anomalousNodes = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				if (!sCurrentLine.equals("null")) {
					Integer lineInt = Integer.parseInt(sCurrentLine);
					anomalousNodes = ArrayUtils.add(anomalousNodes, lineInt);
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
		return anomalousNodes;
	}

	/**
	 * get all the nodes in the fileName, it will return Integer array
	 *
	 * @return all the nodes in the fileName
	 */
	public static double[] getCounts(String fileName) {
		BufferedReader br = null;
		double[] anomalousNodes = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				if (!sCurrentLine.equals("null")) {
					String[] lines = sCurrentLine.split(" ");
					Integer lineInt = Integer.parseInt(lines[1]);
					anomalousNodes = ArrayUtils.add(anomalousNodes, lineInt);
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
		return anomalousNodes;
	}

	/**
	 * return the nodes in the file with int array
	 *
	 * @return all the nodes in the fileName
	 */
	public static int[] getIntFromFile(String fileName) {
		BufferedReader br = null;
		int[] anomalousNodes = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				if (!sCurrentLine.equals("null")) {
					Integer lineInt = Integer.parseInt(sCurrentLine);
					anomalousNodes = ArrayUtils.add(anomalousNodes, lineInt);
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
		return anomalousNodes;
	}

	/**
	 * return the nodes in the file with int array
	 *
	 * @return all the nodes in the fileName
	 */
	public static int[] getIntFromFile(File file, boolean flag) {
		BufferedReader br = null;
		int[] anomalousNodes = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(file));
			while ((sCurrentLine = br.readLine()) != null) {
				if (!sCurrentLine.equals("null")) {
					Integer lineInt = Integer.parseInt(sCurrentLine);
					anomalousNodes = ArrayUtils.add(anomalousNodes, lineInt);
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
		return anomalousNodes;
	}

	/**
	 * return the nodes in the file with int array
	 *
	 * @return all the nodes in the fileName
	 */
	public static long[] getIntFromFile(File file) {
		BufferedReader br = null;
		long[] anomalousNodes = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(file));
			while ((sCurrentLine = br.readLine()) != null) {
				if (!sCurrentLine.equals("null")) {
					Long lineInt = Long.parseLong(sCurrentLine);
					anomalousNodes = ArrayUtils.add(anomalousNodes, lineInt);
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
		return anomalousNodes;
	}

	/**
	 * return the nodes in the file with int array
	 *
	 * @return all the nodes in the fileName
	 */
	public static int[] getCiCntFromFile(String fileName) {
		BufferedReader br = null;
		int[] anomalousNodes = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				if (!sCurrentLine.equals("null")) {
					String[] tmp = sCurrentLine.split(" ");
					Integer lineInt = Integer.parseInt(tmp[1]);
					anomalousNodes = ArrayUtils.add(anomalousNodes, lineInt);
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
		return anomalousNodes;
	}

	public static HashMap<Integer, Integer> getCitCntsFromFile(String fileName) {
		HashMap<Integer, Integer> hashMap = new HashMap<Integer, Integer>();
		BufferedReader br = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				String[] str = sCurrentLine.split(" ");
				int node = Integer.parseInt(str[0]);
				Integer pValue = Integer.parseInt(str[2]);
				hashMap.put(node, pValue);
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
		return hashMap;
	}

	/**
	 * return the nodes in the file with int array
	 *
	 * @return all the nodes in the fileName
	 */
	public static int[] getWaterDataNodes(String fileName) {
		BufferedReader br = null;
		int[] anomalousNodes = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				if (!sCurrentLine.equals("null") && !sCurrentLine.isEmpty()) {
					String[] lineInfo = sCurrentLine.split(" ");
					if (lineInfo[0].trim().equals("Junc") || lineInfo[0].trim().equals("Tank")
							|| lineInfo[0].trim().equals("RESERVOIR")) {
						String[] nodeInfo = lineInfo[1].split("-");
						Integer lineInt = Integer.parseInt(nodeInfo[1]);
						double chemicalValue = Double.parseDouble(lineInfo[lineInfo.length - 1].trim());

						if (chemicalValue >= 0.8) {
							anomalousNodes = ArrayUtils.add(anomalousNodes, lineInt);
						}
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
		return anomalousNodes;
	}

	/**
	 * return the nodes in the file with int array
	 *
	 * @return all the nodes in the fileName
	 */
	public static double[] getPValueFromEcuador(String fileName) {
		BufferedReader br = null;
		double[] PValue = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				String[] lineInfo = sCurrentLine.split(" ");
				PValue = ArrayUtils.add(PValue, Double.parseDouble(lineInfo[1].trim()));
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
		return PValue;
	}

	/**
	 * read all the double value from the fileName, it will return Double array
	 *
	 * @return read all the double value from the fileName
	 */
	public static Double[] getDoubleFromFile(String fileName) {
		BufferedReader br = null;
		Double[] doubleValues = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				Double lineInt = Double.parseDouble(sCurrentLine);
				doubleValues = ArrayUtils.add(doubleValues, lineInt);
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
		return doubleValues;
	}

	/**
	 * read all the double value from the fileName, it will return double array
	 */
	public static double[] getdoubleFromFile(String fileName) {
		BufferedReader br = null;
		double[] doubleValues = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				double lineInt = Double.parseDouble(sCurrentLine);
				doubleValues = ArrayUtils.add(doubleValues, lineInt);
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
		return doubleValues;
	}

	public static HashMap<Integer, Double> getPValueFromFile(String fileName) {
		HashMap<Integer, Double> hashMap = new HashMap<Integer, Double>();
		BufferedReader br = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				String[] str = sCurrentLine.split(" ");
				int node = Integer.parseInt(str[0]);
				double pValue = Double.parseDouble(str[1]);
				hashMap.put(node, pValue);
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
		return hashMap;
	}

	/**
	 * the file which contains all the edges, each edge shows two times we use
	 * the "\t" or " " to split the edge end points
	 */
	public static ArrayList<ArrayList<Integer>> getGraphFromFile(String fileName, String splitter) {
		BufferedReader br = null;
		ArrayList<ArrayList<Integer>> graphList = new ArrayList<ArrayList<Integer>>();
		int count = 0;
		int flag = 0;

		try {
			String sCurrentLine;

			br = new BufferedReader(new FileReader(fileName));
			ArrayList<Integer> arrayList = new ArrayList<Integer>();
			while ((sCurrentLine = br.readLine()) != null) {
				String[] str = sCurrentLine.split(splitter);
				Integer[] lineInt = new Integer[2];
				lineInt[0] = Integer.parseInt(str[0]);
				lineInt[1] = Integer.parseInt(str[1]);
				if (count == 0) {
					flag = lineInt[0];
					arrayList.add(lineInt[1]);
				} else {
					if (flag == lineInt[0]) {
						arrayList.add(lineInt[1]);
					} else {
						graphList.add(arrayList);
						arrayList = new ArrayList<Integer>();
						flag = lineInt[0];
						arrayList.add(lineInt[1]);
					}
				}
				count++;
			}
			graphList.add(arrayList); // add the last line .
			if (!graphList.isEmpty()) {
				System.out.println("\n======================================================="
						+ "\ntotal number of nodes in graph is : " + graphList.size()
						+ "\ntotal number of the edges in graph is : " + (int) (count / 2)
						+ "\n=======================================================\n");
			} else {
				System.out.println("there does not exist any information in the file!!");
				System.exit(0);
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
		return graphList;
	}

	public static ArrayList<Integer[]> getGraphEdgeFromFile(String fileName, String splitter) {
		ArrayList<Integer[]> arr = new ArrayList<Integer[]>();
		BufferedReader br = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				String[] str = sCurrentLine.split(splitter);
				Integer[] lineInt = new Integer[2];
				lineInt[0] = Integer.parseInt(str[0]);
				lineInt[1] = Integer.parseInt(str[1]);
				arr.add(lineInt);
			}
			if (!arr.isEmpty()) {
				System.out.println("total number of the graph edges is : " + arr.size());
			} else {
				System.out.println("there does not exist any information in the file!!");
				System.exit(0);
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
		return arr;
	}

	public static ArrayList<Edge> getGraphEdgeFromFile(String fileName) {
		ArrayList<Edge> arr = new ArrayList<Edge>();
		BufferedReader br = null;
		try {
			String sCurrentLine;
			int count = 0;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				String[] str = sCurrentLine.split(" ");
				Integer[] lineInt = new Integer[2];
				lineInt[0] = Integer.parseInt(str[0]);
				lineInt[1] = Integer.parseInt(str[1]);
				arr.add(new Edge(lineInt[0], lineInt[1], count++, Double.parseDouble(str[2])));
			}
			if (!arr.isEmpty()) {
				System.out.println("total number of the graph edges is : " + arr.size());
			} else {
				System.out.println("there does not exist any information in the file!!");
				System.exit(0);
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
		return arr;
	}

	public static HashMap<String, Double> loadRealDataGraphFromEmbers(String fileName, String splitter) {
		HashMap<String, Double> graph = new HashMap<String, Double>();
		BufferedReader br = null;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			while ((sCurrentLine = br.readLine()) != null) {
				String[] str = sCurrentLine.split(splitter);
				String edgeName = str[0];
				Double distance = Double.parseDouble(str[1]);
				graph.put(edgeName, distance);
			}
			if (!graph.isEmpty()) {
				System.out.println("total number of the graph edges is : " + graph.size());
			} else {
				System.out.println("there does not exist any information in the file!!");
				System.exit(0);
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
		return graph;
	}

	/**
	 * save Integer nodes to the fileName
	 */
	public static void saveIntegerNodes2File(Integer[] nodes, String fileName) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName);
		String str = "";
		if (nodes != null) {
			for (int i = 0; i < nodes.length; i++) {
				str = str + nodes[i];
				str = str + "\n";
			}
		} else {
			str = "null";
		}
		fw.write(str);
		fw.close();
	}

	/**
	 * save int nodes to the fileName
	 */
	public static void saveIntNodes2File(int[] nodes, String fileName) {
		FileWriter fw = null;
		try {
			fw = new FileWriter(fileName);
			String str = "";
			if (nodes != null) {
				for (int i = 0; i < nodes.length; i++) {
					str = str + nodes[i];
					str = str + "\n";
				}
			} else {
				str = "null";
			}
			fw.write(str);
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * save int nodes to the fileName
	 */
	public static void saveIntNodes2File(long[] nodes, String fileName) {
		FileWriter fw = null;
		try {
			fw = new FileWriter(fileName);
			String str = "";
			if (nodes != null) {
				for (int i = 0; i < nodes.length; i++) {
					str = str + nodes[i];
					str = str + "\n";
				}
			} else {
				str = "null";
			}
			fw.write(str);
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * save int nodes to the fileName
	 */
	public static void saveIntNodes2File(int[] nodes, String fileName, double score) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName);
		String date = fileName.split("\\.")[0];
		String str = "";
		if (nodes != null) {
			for (int i = 0; i < nodes.length; i++) {
				str = str + nodes[i] + " " + String.valueOf(score) + " " + date + "\n";
				str = str + "\n";
			}
		} else {
			str = "null";
		}
		fw.write(str);
		fw.close();
	}

	/**
	 * save graph edge to the fileName
	 */
	public static void saveEdges2File(ArrayList<Integer[]> edges, String fileName, double totalScore, int rootID,
			boolean flag) throws IOException {
		FileWriter fw = null;
		fileName = "realDataResult/ModModPCST/" + fileName;
		fw = new FileWriter(fileName);
		String str = "";
		str = str + totalScore + "\n";
		if (edges != null) {
			for (Integer[] e : edges) {
				str = str + e[0] + "\t" + e[1];
				str = str + "\n";
			}
		} else {
			str = str + rootID;
		}
		fw.write(str);
		fw.close();
	}

	/**
	 * save graph edge to the fileName
	 */
	public static void saveEdges2File(ArrayList<Integer[]> edges, String fileName, boolean flag) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName);
		String str = "";
		if (edges != null) {
			for (Integer[] e : edges) {
				str = str + e[0] + " " + e[1];
				str = str + "\n";
			}
		} else {
			str = null;
		}
		fw.write(str);
		fw.close();
	}

	public static void saveEdges2File(ArrayList<Edge> edges, String fileName) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName);
		String str = "";
		if (edges != null) {
			for (Edge e : edges) {
				str = str + e.i + "\t" + e.j;
				str = str + "\n";
			}
		} else {
			str = "null";
		}
		fw.write(str);
		fw.close();
	}

	public static void savePValue2File(double[] pValue, String fileName) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName);
		String str = "";
		if (pValue != null) {
			for (int i = 0; i < pValue.length; i++) {
				str = str + pValue[i];
				str = str + "\n";
			}
		} else {
			str = "null";
		}
		fw.write(str);
		fw.close();
	}

	public static void savePValue2File(double[] pValue, String fileName, boolean flag) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName);
		String str = "";
		if (pValue != null) {
			for (int i = 0; i < pValue.length; i++) {
				str = str + i + "\t" + pValue[i];
				str = str + "\n";
			}
		} else {
			str = "null";
		}
		fw.write(str);
		fw.close();
	}

	/**
	 * create a new file and save the doubleArray to the fileName
	 */
	public static void saveDoubleArray2File(Double[] doubleArray, String fileName) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName);
		String str = "";
		for (int i = 0; i < doubleArray.length; i++) {
			str = str + doubleArray[i];
			str = str + "\n";
		}
		fw.write(str);
		fw.close();
	}

	/**
	 * save doubleArray to fileName, if the file exists, we just append the
	 * double array to the file
	 */
	public static void saveDoubleArray2FileAppend(Double[] doubleArray, String fileName) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName, true);
		String str = "";
		for (int i = 0; i < doubleArray.length; i++) {
			str = str + doubleArray[i];
			str = str + "\n";
		}
		fw.write(str);
		fw.close();
	}

	/**
	 * save doubleArray to fileName, if the file exists, we just append this
	 * array to the fileName
	 */
	public static void savedoubleArray2FileAppend(double[] doubleArray, String fileName) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName, true);
		String str = "";
		for (int i = 0; i < doubleArray.length; i++) {
			str = str + doubleArray[i];
			str = str + "\n";
		}
		fw.write(str);
		fw.close();
	}

	public static void saveData(String fileName, double totalTime, int[] resultSubset) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName);
		String str1 = "totalTime is:" + totalTime + "\n";
		String str2 = "resultSubset: " + "\n";
		fw.write(str1);
		fw.write(str2);
		String str = "";
		for (int i = 0; i < resultSubset.length; i++) {
			str = str + i + " ";
			if (i % 30 == 0 && i != 0)
				str = str + "\n";
		}
		fw.write(str);
		str = "\nlength of result subset is:" + resultSubset.length + "\n";
		fw.write(str);
		fw.close();
	}

	public static void saveData(String fileName, double totalTime, Integer[] resultSubset) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName);
		String str1 = "totalTime is:" + totalTime + "\n";
		String str2 = "resultSubset: " + "\n";
		fw.write(str1);
		fw.write(str2);
		String str = "";
		for (int i = 0; i < resultSubset.length; i++) {
			str = str + i + " ";
			if (i % 30 == 0 && i != 0)
				str = str + "\n";
		}
		fw.write(str);
		str = "\nlength of result subset is:" + resultSubset.length + "\n";
		fw.write(str);
		fw.close();
	}

	/**
	 * save double array ; if there is no file, it creates this fileName,
	 * otherwise it will override this file
	 */
	public static void savedoubleArray2File(double[] abnormalNodes, String fileName) throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName, false);
		String str = "";
		for (int i = 0; i < abnormalNodes.length; i++) {
			str = str + abnormalNodes[i];
			str = str + "\n";
		}
		fw.write(str);
		fw.close();
	}

	/**
	 * @param flag
	 *            can control append the array to the fileName or not
	 */
	public static void savedoubleArray2FileAppend(double[] doubleArray, String fileName, boolean flag)
			throws IOException {
		FileWriter fw = null;
		fw = new FileWriter(fileName, flag);
		String str = "";
		for (int i = 0; i < doubleArray.length; i++) {
			str = str + doubleArray[i];
			str = str + "\n";
		}
		fw.write(str);
		fw.close();
	}

	public static double[] generatePValue(int size, boolean flag, double alphaMax) {
		if (flag == false) {
			if (size > 0) {
				double[] PValue = new double[size];
				return PValue;
			} else {
				return null;
			}
		} else {
			double[] PValue = new double[100];
			int[] anomalousNodes = new int[] { 11, 12, 13, 21, 23, 31, 32, 33, 16, 17, 18, 26, 27, 28, 36, 37, 38, 54,
					55, 64, 65, 74, 75, 84, 85, 99 };
			for (int i = 0; i < PValue.length; i++) {
				if (ArrayUtils.contains(anomalousNodes, i)) {
					Random random = new Random();
					double value = random.nextDouble();
					while (value > alphaMax) {
						value = random.nextDouble();
					}
					PValue[i] = value;
				} else {
					Random random = new Random();
					double value = random.nextDouble();
					while (value < alphaMax) {
						value = random.nextDouble();
					}
					PValue[i] = value;
				}
			}
			return PValue;
		}
	}

	/**
	 * simulation of the AdjMatrix for grid data in this simulation graph, each
	 * of the node has 4 neighbors
	 *
	 * @param N
	 *            we have N*N number of nodes need to simulate
	 */
	public static boolean[][] generateGridGraph(int N) {
		boolean[][] graph = new boolean[N * N][N * N];
		int[][] Nodes = new int[N][N];
		int count = 0;
		for (int i = 0; i < Nodes.length; i++)
			for (int j = 0; j < Nodes.length; j++) {
				Nodes[i][j] = count++;
			}
		count = 0;
		int i = 1;
		for (int k = 1; k <= N * N; k++) {
			count++;
			if (count > N) {
				count = 1;
				i++;
			}
			int j = count;
			if ((i - 1) < 1)
				graph[k - 1][Nodes[N - 1][j - 1]] = true;
			else
				graph[k - 1][Nodes[i - 2][j - 1]] = true;
			if ((i + 1) > N)
				graph[k - 1][Nodes[0][j - 1]] = true;
			else
				graph[k - 1][Nodes[i][j - 1]] = true;
			if ((j - 1) < 1)
				graph[k - 1][Nodes[i - 1][N - 1]] = true;
			else
				graph[k - 1][Nodes[i - 1][j - 2]] = true;
			if ((j + 1) > N)
				graph[k - 1][Nodes[i - 1][0]] = true;
			else
				graph[k - 1][Nodes[i - 1][j]] = true;
		}
		return graph;
	}

	public static double[][] generateGraph(int N, double weight) {
		double[][] graph = new double[N * N][N * N];
		for (int k = 0; k < graph.length; k++) {
			for (int j = 0; j < graph.length; j++) {
				graph[k][j] = -1;
			}
		}
		int[][] Nodes = new int[N][N];
		int count = 0;
		for (int i = 0; i < Nodes.length; i++)
			for (int j = 0; j < Nodes.length; j++) {
				Nodes[i][j] = count++;
			}
		count = 0;
		int i = 1;
		for (int k = 1; k <= N * N; k++) {
			count++;
			if (count > N) {
				count = 1;
				i++;
			}
			int j = count;
			if ((i - 1) < 1)
				graph[k - 1][Nodes[N - 1][j - 1]] = weight;
			else
				graph[k - 1][Nodes[i - 2][j - 1]] = weight;
			if ((i + 1) > N)
				graph[k - 1][Nodes[0][j - 1]] = weight;
			else
				graph[k - 1][Nodes[i][j - 1]] = weight;
			if ((j - 1) < 1)
				graph[k - 1][Nodes[i - 1][N - 1]] = weight;
			else
				graph[k - 1][Nodes[i - 1][j - 2]] = weight;
			if ((j + 1) > N)
				graph[k - 1][Nodes[i - 1][0]] = weight;
			else
				graph[k - 1][Nodes[i - 1][j]] = weight;
		}
		return graph;
	}

	public static double[][] generateGraph2(int N) {
		double[][] graph = new double[N * N][N * N];
		for (int k = 0; k < graph.length; k++) {
			for (int j = 0; j < graph.length; j++) {
				graph[k][j] = -1;
			}
		}
		int[][] Nodes = new int[N][N];
		int count = 0;
		for (int i = 0; i < Nodes.length; i++)
			for (int j = 0; j < Nodes.length; j++) {
				Nodes[i][j] = count++;
			}
		count = 0;
		int i = 1;
		for (int k = 1; k <= N * N; k++) {
			count++;
			if (count > N) {
				count = 1;
				i++;
			}
			int j = count;
			if ((i - 1) < 1)
				graph[k - 1][Nodes[N - 1][j - 1]] = 1;
			else
				graph[k - 1][Nodes[i - 2][j - 1]] = 1;
			if ((i + 1) > N)
				graph[k - 1][Nodes[0][j - 1]] = 1;
			else
				graph[k - 1][Nodes[i][j - 1]] = 1;
			if ((j - 1) < 1)
				graph[k - 1][Nodes[i - 1][N - 1]] = 1;
			else
				graph[k - 1][Nodes[i - 1][j - 2]] = 1;
			if ((j + 1) > N)
				graph[k - 1][Nodes[i - 1][0]] = 1;
			else
				graph[k - 1][Nodes[i - 1][j]] = 1;
		}
		return graph;
	}

	/**
	 * according to the graph, we generate subgraph using the random walk. for
	 * example, we want to select 20% of nodes as event, we need to set percent
	 * to 0.2
	 */
	public static int[] generateSubGraph(double[][] graph, double percent) {
		int[] result = null;
		if (percent > 0.5) {
			System.out.println("the percentage is too large.");
			System.exit(0);
		}
		int size = (int) percent * graph.length;
		Random random = new Random();
		int startNode = random.nextInt(graph.length);
		boolean flag = true;
		while (flag) {
			result = ArrayUtils.add(result, startNode);
			int[] adjNodes = null;
			for (int i = 0; i < graph.length; i++) {
				if (graph[startNode][i] >= 0) {
					adjNodes = ArrayUtils.add(adjNodes, i);
				}
			}
			startNode = adjNodes[random.nextInt(adjNodes.length)];
			if (result.length >= size) {
				break;
			}
		}
		return result;
	}

	public static int[] generateSubGraph(double[][] graph, int numSubGraph) {
		int[] result = null;
		for (int i = 0; i < numSubGraph; i++) {
			result = ArrayUtils.addAll(result, Utils.generateSubGraph(graph, 0.1));
		}
		HashSet<Integer> hashSet = new HashSet<Integer>();
		for (Integer i : result) {
			hashSet.add(i);
		}
		result = null;
		for (Integer i : hashSet) {
			result = ArrayUtils.add(result, i);
		}
		return result;
	}

	public static boolean[][] generateUnWeightedGraph(ArrayList<ArrayList<Integer>> graphList) {
		boolean[][] graph = new boolean[graphList.size()][graphList.size()];
		for (int i = 0; i < graph.length; i++)
			for (int j = 0; j < graph.length; j++)
				graph[i][j] = false;
		int count = 0;
		for (int i = 0; i < graphList.size(); i++) {
			for (int j = 0; j < graphList.get(i).size(); j++) {
				graph[i][graphList.get(i).get(j)] = true;
				count++;
			}
		}
		System.out.println("total number of the edges is: " + count / 2);
		return graph;
	}

	/**
	 * we model our graph's edge weight with different situations. the graph
	 * constructs like this: if graph[i][j] larger or equal to 0 stands for
	 * there exists an edge. else there is no edge in node i and node j.
	 */
	public static double[][] generateWeightedGraph(boolean[][] graph, double[] pi, int[] abnormalNodes) {
		double[][] weightedGraph = new double[graph.length][graph.length];
		for (int i = 0; i < weightedGraph.length; i++)
			for (int j = 0; j < weightedGraph.length; j++)
				weightedGraph[i][j] = -1;
		for (int i = 0; i < graph.length; i++) {
			for (int j = i + 1; j < graph.length; j++) {
				// if there is an edge, we will process this edge
				if (graph[i][j] == true) {
					// if two vertices pi value larger than 0, we let weight to
					// 0
					if (pi[i] > 0 && pi[j] > 0) {
						weightedGraph[i][j] = 0;
						weightedGraph[j][i] = 0;
						// if one vertex pi value larger than 0 and another
						// vertex pi value less than 0, we let process as
						// follows
					} else if (pi[i] > 0 && pi[j] < 0) {
						int count = 0;
						for (int k : abnormalNodes) {
							if (graph[k][j] == true) {
								count++;
							}
						}
						if (count != 0) {
							double weightValue = -pi[j] / count;
							weightedGraph[i][j] = weightValue;
							weightedGraph[j][i] = weightValue;
						} else {
							System.out.println("something was wrong");
							System.exit(0);
						}
						// if one vertex pi value larger than 0 and another
						// vertex pi value less than 0, we let process as
						// follows
					} else if (pi[i] < 0 && pi[j] > 0) {
						int count = 0;
						for (int k : abnormalNodes) {
							if (graph[k][i] == true) {
								count++;
							}
						}
						if (count != 0) {
							double weightValue = -pi[i] / count;
							weightedGraph[i][j] = weightValue;
							weightedGraph[j][i] = weightValue;
						} else {
							System.out.println("something was wrong");
							System.exit(0);
						}
						// if two vertices are less than 0, the weight will be
						// the average of two vertices value.
					} else {
						if (pi[i] == 0 || pi[j] == 0) {
							System.out.println("there exists a bug");
							System.exit(0);
						}
						weightedGraph[i][j] = -(pi[i] + pi[j]) / 2;
						weightedGraph[j][i] = -(pi[i] + pi[j]) / 2;
					}
				}
			}
		}
		System.out.println("generated weighted graph finish");
		return weightedGraph;
	}

	/**
	 * if a == b return true; else return false
	 */
	public static boolean compareTwoArray(Integer[] a, Integer[] b) {
		boolean flag = true;
		if (a == null && b == null) {
			flag = false;
		} else if (a != null && b != null) {
			if (a.length != b.length) {
				flag = false;
			} else {
				Arrays.sort(a);
				Arrays.sort(b);
				for (int i = 0; i < a.length; i++) {
					if (a[i].intValue() != b[i].intValue()) {
						flag = false;
						break;
					}
				}
			}
		} else {
			flag = false;
		}
		return flag;
	}

	/**
	 * if a == b return true; else return false
	 */
	public static boolean compareTwoArray(int[] a, int[] b) {
		boolean flag = true;
		if (a == null && b == null) {
			flag = false;
		} else if (a != null && b != null) {
			if (a.length != b.length) {
				flag = false;
			} else {
				Arrays.sort(a);
				Arrays.sort(b);
				for (int i = 0; i < a.length; i++) {
					if (a[i] != b[i]) {
						flag = false;
						break;
					}
				}
			}
		} else {
			flag = false;
		}
		return flag;
	}

	/**
	 * convert Integer to int array
	 *
	 * @return the int array
	 */
	public static int[] integerToInt(Integer[] array) {
		int[] resultArray = null;
		if (array != null) {
			for (int i : array) {
				resultArray = ArrayUtils.add(resultArray, i);
			}
		} else {
			return null;
		}

		return resultArray;
	}

	public static Integer[] intToInteger(int[] array) {
		if (array == null) {
			return null;
		} else {
			Integer[] newArray = new Integer[array.length];
			for (int i = 0; i < newArray.length; i++) {
				newArray[i] = array[i];
			}
			return newArray;
		}
	}

	public static void generateSTPFile(String name, ArrayList<Edge> edges, double[] PValue) {
		DecimalFormat decimalFormat = new DecimalFormat("#.######");
		String strHead = "33D32945 STP File, STP Format Version 1.0\n" + "\n" + "SECTION Comment\n" + "Name \"" + name
				+ "\"\n" + "Creator \"Baojian \"" + "Problem \"Prize-collecting Steiner Tree\"\n"
				+ "Remark \"Automatically converted from MWCS\"\n" + "END\n" + "\n" + "SECTION Graph";
		strHead = "Nodes " + PValue.length + "\n";
		strHead = "Edges " + edges.size() + "\n";

		if (edges != null) {
			for (Edge e : edges) {
				strHead = strHead + "E " + e.i + " " + e.j + " " + decimalFormat.format(e.cost) + "\n";
			}
		}

		strHead = strHead + "END\n\n" + "SECTION Terminals\n";
		if (PValue != null) {
			strHead = strHead + "Terminals " + PValue.length + "\n";
			for (int i = 0; i < PValue.length; i++) {
				strHead = strHead + "TP " + i + " " + decimalFormat.format(PValue[i]);
			}
		}
		System.out.println(strHead);
	}

	public static void main(String args[]) {
		generateSTPFile("gridData", null, null);
	}

	public static double[] getChemical(String fileName) {
		BufferedReader br = null;
		double[] chemical = new double[12527];
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(fileName));
			int count = 0;
			while ((sCurrentLine = br.readLine()) != null) {
				if (!sCurrentLine.equals("null") && !sCurrentLine.isEmpty()) {
					String[] lineInfo = sCurrentLine.split(" ");
					if (lineInfo[0].trim().equals("Junc") || lineInfo[0].trim().equals("Tank")
							|| lineInfo[0].trim().equals("RESERVOIR")) {
						chemical[count++] = Double.parseDouble(lineInfo[lineInfo.length - 1].trim());
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
		return chemical;
	}

	public static HashMap<String, ArrayList<FuncType>> checkAlreadyDone(String resultFileName) {
		try {
			if (!new File(resultFileName).exists()) {
				return null;
			}
			HashMap<String, ArrayList<FuncType>> done = new HashMap<String, ArrayList<FuncType>>();
			BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(resultFileName)));
			String line = null;
			while ((line = bufferedReader.readLine()) != null) {
				line = line.trim();
				String[] items = line.split(";");
				String fileName = items[0].trim();
				FuncType funcID = FuncType.valueOf(items[2].trim());
				if (done.containsKey(fileName)) {
					done.get(fileName).add(funcID);
				} else {
					ArrayList<FuncType> arr = new ArrayList<FuncType>();
					arr.add(funcID);
					done.put(fileName, arr);
				}
			}
			bufferedReader.close();
			return done;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}

	public static int[] toArray(HashSet<Integer> set) {
		int[] result = null;
		for (int i : set) {
			result = ArrayUtils.add(result, i);
		}
		return result;
	}

}
