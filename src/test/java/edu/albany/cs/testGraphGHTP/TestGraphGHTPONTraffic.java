package edu.albany.cs.testGraphGHTP;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.lang3.ArrayUtils;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.Stat;
import edu.albany.cs.graphGHTP.GraphGHTP;
import edu.albany.cs.scoreFuncs.FuncType;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.ScoreFuncFactory;

public class TestGraphGHTPONTraffic {
	private final int numOfThreads;

	private final String trafficDataFolder = "data/Traffic/";
	private final String resultFolder = "./output/Traffic/";
	private final String resultFileName = resultFolder + "graph_GHTP_Traffic_Result.txt";
	private String[] testingDates = new String[] { "2014-03-01", "2014-03-02", "2014-03-03", "2014-03-04", "2014-03-05",
			"2014-03-06", "2014-03-07", "2014-03-08", "2014-03-09", "2014-03-10", "2014-03-11", "2014-03-12",
			"2014-03-13", "2014-03-14", "2014-03-15", "2014-03-16", "2014-03-17", "2014-03-18", "2014-03-19",
			"2014-03-20", "2014-03-21", "2014-03-22", "2014-03-23", "2014-03-24", "2014-03-25", "2014-03-26",
			"2014-03-27", "2014-03-28", "2014-03-29", "2014-03-30", "2014-03-31", };

	private int verboseLevel = 0;

	public TestGraphGHTPONTraffic(int numOfThreads) {
		testingDates = new String[] { "2014-03-01", "2014-03-02", "2014-03-03", "2014-03-04", "2014-03-05",
				"2014-03-06", };
		this.numOfThreads = numOfThreads;
		run();
	}

	private void run() {

		if (!new File(resultFolder).isDirectory()) {
			new File(resultFolder).mkdir();
		}

		String[] subFolders = new String[testingDates.length];
		for (int i = 0; i < testingDates.length; i++) {
			subFolders[i] = trafficDataFolder + testingDates[i];
		}
		for (int i = 0; i < testingDates.length; i++) {

			File[] allFilesInSubFolder = new File(subFolders[i]).listFiles();
			final CountDownLatch latch = new CountDownLatch(allFilesInSubFolder.length);
			ExecutorService pool = Executors.newFixedThreadPool(numOfThreads);
			for (final File apdmFile : allFilesInSubFolder) {

				final APDMInputFormat apdm = new APDMInputFormat(apdmFile);
				final int graphSize = apdm.data.numNodes;
				final ArrayList<Integer[]> edges = apdm.data.intEdges;
				final ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;

				pool.execute(new Thread() {
					public void run() {

						System.out.println("processing : " + apdmFile.getName() + " funcID: " + FuncType.EMS);
						long startTime = System.nanoTime();
						Function func = null;
						double[] b = new double[apdm.data.numNodes];
						double[] PValue = apdm.data.PValue;
						int[] trueSubGraph = apdm.data.trueSubGraphNodes;

						double[] logC = new double[apdm.data.numNodes];
						for (int i = 0; i < apdm.data.numNodes; i++) {
							logC[i] = -Math.log(PValue[i] / 0.15D);
						}
						double mean = new Stat(logC).mean();
						double std = new Stat(logC).std();
						for (int i = 0; i < apdm.data.numNodes; i++) {
							if(std != 0.0D){
								logC[i] = (logC[i] - mean) / std;	
							}
						}
						Arrays.fill(b, 1.0D);
						func = ScoreFuncFactory.getFunc(FuncType.EMS, b, logC);
						double bestFuncValue = -Double.MAX_VALUE;
						double[] bestFuncs = null;
						GraphGHTP bestGraphGHTP = null;
						for (int s = 50; s <= 1000; s += 50) {
							int g = 1;
							double B = s - 1 + 0.0D;
							boolean isInitialSingle = false;
							GraphGHTP graphGHTP = new GraphGHTP(graphSize, edges, edgeCosts, logC, s, g, B,
									isInitialSingle, trueSubGraph, func, true, 1.0D);
							if (bestFuncValue < graphGHTP.funcValueTail) {
								bestFuncValue = graphGHTP.funcValueTail;
								bestGraphGHTP = graphGHTP;
								bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
							}
							if (verboseLevel == 0) {
								System.out.println("processing s: " + s + " ; score: " + graphGHTP.funcValueTail);
							}
						}
						try {
							double runTime = (System.nanoTime() - startTime) / 1e9D;
							double funcVal = bestGraphGHTP.funcValueTail;
							int numOfNodes = bestGraphGHTP.resultNodesTail.size();
							int numOfNodesP = 0;
							int numOfNOdesLargerThan0 = 0;
							for (int node : bestGraphGHTP.resultNodesTail) {
								if (PValue[node] <= 0.15D) {
									numOfNodesP++;
								}
								if (logC[node] > 0.0D) {
									numOfNOdesLargerThan0++;
								}
							}
							FileWriter fileWriter = new FileWriter(new File(resultFileName), true);
							fileWriter.write(apdmFile.getName() + "," + FuncType.EMS + "," + funcVal + "," + runTime
									+ ",0.0,0.0,0.0,0.0 0.0," + numOfNodes + "," + numOfNodesP + ","
									+ numOfNOdesLargerThan0 + "\n");
							fileWriter.close();
						} catch (IOException e) {
							e.printStackTrace();
						}
						latch.countDown();
					}

				});
			}
			try {
				latch.await();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			pool.shutdown();
		}
	}

	public static void main(String args[]) {

		if (args == null || args.length == 0) {
			int numOfThreads = 3;
			new TestGraphGHTPONTraffic(numOfThreads);
		} else {
			int numOfThreads = Integer.parseInt(args[0]);
			new TestGraphGHTPONTraffic(numOfThreads);
		}
	}

}
