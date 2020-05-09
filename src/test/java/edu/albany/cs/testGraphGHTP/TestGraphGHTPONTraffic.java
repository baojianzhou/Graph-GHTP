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
import edu.albany.cs.base.Constants;
import edu.albany.cs.base.Stat;
import edu.albany.cs.graphGHTP.GraphGHTP;
import edu.albany.cs.scoreFuncs.FuncType;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.ScoreFuncFactory;
import org.junit.Test;

public class TestGraphGHTPONTraffic {

	private final String resultFileName = Constants.TrafficOutputFolder + "graph_GHTP_Traffic_Result.txt";

	private int verboseLevel = 0;

	private void run(int numOfThreads, String[] testingDates) {

		String[] subFolders = new String[testingDates.length];
		for (int i = 0; i < testingDates.length; i++) {
			subFolders[i] = Constants.TrafficDataFolder + testingDates[i];
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

	@Test
	public void runTest() {
		
		Constants.intializeProject();
		int numOfThreads = 3;
		run(numOfThreads, new String[] { "2014-03-01", "2014-03-02", "2014-03-03", "2014-03-04", "2014-03-05",
				"2014-03-06", });
	}

}
