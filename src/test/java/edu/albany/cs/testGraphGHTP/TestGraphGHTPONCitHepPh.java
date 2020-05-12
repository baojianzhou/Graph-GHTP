package edu.albany.cs.testGraphGHTP;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.Constants;
import edu.albany.cs.graphGHTP.GraphGHTP;
import edu.albany.cs.scoreFuncs.FuncType;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.ScoreFuncFactory;
import org.junit.Test;

public class TestGraphGHTPONCitHepPh {

	private final String resultFileName = Constants.CitHepPhOutputFolder + "graph_GHTP_CitHepPh_Result_Updated.txt";


	private void run(int numOfThreads) {

		ExecutorService pool = Executors.newFixedThreadPool(numOfThreads);

		for (final File apdmFile : new File(Constants.CitHepPhDataFolder).listFiles()) {

			final APDMInputFormat apdm = new APDMInputFormat(apdmFile);
			final int graphSize = apdm.data.numNodes;
			final ArrayList<Integer[]> edges = apdm.data.intEdges;
			final ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;

			for (final FuncType funcID : Constants.funcs) {

				pool.execute(new Thread() {
					public void run() {
						System.out.println("processing : " + apdmFile.getName() + " funcID: " + funcID);
						long startTime = System.nanoTime();
						Function func = null;
						double[] b = new double[apdm.data.numNodes];
						double[] c = apdm.data.counts;
						int[] trueSubGraph = apdm.data.trueSubGraphNodes;
						switch (funcID) {
						case Kulldorff:
							c = apdm.data.counts;
							double maximalCount = StatUtils.max(apdm.data.counts);
							double meanCount = StatUtils.mean(apdm.data.counts);
							System.out.println("maximalValue: " + maximalCount);
							Arrays.fill(b, maximalCount);
							Arrays.fill(b, meanCount);
							func = ScoreFuncFactory.getFunc(funcID, b, c);

							double bestFuncValue = -Double.MAX_VALUE;
							double[] bestFuncs = null;
							GraphGHTP bestGraphGHTP = null;

							for (int s = 50; s <= 1000; s += 50) {
								int g = 1;
								double B = s - 1 + 0.0D;
								GraphGHTP graphGHTP = new GraphGHTP(graphSize, edges, edgeCosts, apdm.data.counts, s, g,
										B, false, trueSubGraph, func, true, 0.1D);
								int resultSize = graphGHTP.resultNodesTail.size();
								System.out.println("s: " + s + " ; resultSize: " + resultSize + " ; funcValue: "
										+ graphGHTP.funcValueTail);
								if (bestFuncValue < graphGHTP.funcValueTail ) {
									System.out.println("result size: " + resultSize);
									bestFuncValue = graphGHTP.funcValueTail;
									bestGraphGHTP = graphGHTP;
									bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
								}
							}
							try {
								double runTime = (System.nanoTime() - startTime) / 1e9;
								double funcVal = bestGraphGHTP.funcValueTail;
								double pre = 0.0D;
								double rec = 0.0D;
								double fmeasure = 0.0D;
								String fValues = "0.0 0.0";
								int numOfNodes = bestGraphGHTP.resultNodesTail.size();
								FileWriter fileWriter = new FileWriter(new File(resultFileName), true);
								fileWriter.write(apdmFile.getName() + "," + funcID + "," + funcVal + "," + runTime + ","
										+ pre + "," + rec + "," + fmeasure + "," + fValues + "," + numOfNodes + "\n");
								fileWriter.close();
							} catch (IOException e) {
								e.printStackTrace();
							}
							break;
						case EMS:
							c = apdm.data.counts;
							b = apdm.data.averCounts;
							func = ScoreFuncFactory.getFunc(funcID, b, c);

							bestFuncValue = -Double.MAX_VALUE;
							bestFuncs = null;
							bestGraphGHTP = null;
							for (int s = 50; s <= 1000; s += 50) {
								int g = 1;
								double B = s - 1 + 0.0D;
								GraphGHTP graphGHTP = new GraphGHTP(graphSize, edges, edgeCosts, apdm.data.counts, s, g,
										B, false, trueSubGraph, func, true, 1.0D);
								if (bestFuncValue < graphGHTP.funcValueTail ) {
									bestFuncValue = graphGHTP.funcValueTail;
									bestGraphGHTP = graphGHTP;
									bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
								}
							}
							try {
								double runTime = (System.nanoTime() - startTime) / 1e9;
								double funcVal = bestGraphGHTP.funcValueTail;
								double pre = 0.0D;
								double rec = 0.0D;
								double fmeasure = 0.0D;
								String fValues = "";
								for (double fVal : bestGraphGHTP.fValues) {
									fValues += fVal + " ";
								}
								int numOfNodes = bestGraphGHTP.resultNodesTail.size();
								FileWriter fileWriter = new FileWriter(new File(resultFileName), true);
								fileWriter.write(apdmFile.getName() + "," + funcID + "," + funcVal + "," + runTime + ","
										+ pre + "," + rec + "," + fmeasure + "," + fValues + "," + numOfNodes + "\n");
								fileWriter.close();
							} catch (IOException e) {
								e.printStackTrace();
							}
							break;
						case EBP:
							c = apdm.data.counts;
							b = apdm.data.averCounts;
							func = ScoreFuncFactory.getFunc(funcID, b, c);

							bestFuncValue = -Double.MAX_VALUE;
							bestFuncs = null;
							bestGraphGHTP = null;
							for (int s = 50; s <= 1000; s += 50) {
								int g = 1;
								double B = s - 1 + 0.0D;
								GraphGHTP graphGHTP = new GraphGHTP(graphSize, edges, edgeCosts, apdm.data.counts, s, g,
										B, false, trueSubGraph, func, true, 1.0D);
								if (bestFuncValue < graphGHTP.funcValueTail ) {
									bestFuncValue = graphGHTP.funcValueTail;
									bestGraphGHTP = graphGHTP;
									bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
								}
							}
							try {
								double runTime = (System.nanoTime() - startTime) / 1e9;
								double funcVal = bestGraphGHTP.funcValueTail;
								double pre = 0.0D;
								double rec = 0.0D;
								double fmeasure = 0.0D;
								String fValues = "0.0 0.0";
								int numOfNodes = bestGraphGHTP.resultNodesTail.size();
								FileWriter fileWriter = new FileWriter(new File(resultFileName), true);
								fileWriter.write(apdmFile.getName() + "," + funcID + "," + funcVal + "," + runTime + ","
										+ pre + "," + rec + "," + fmeasure + "," + fValues + "," + numOfNodes + "\n");
								fileWriter.close();
							} catch (IOException e) {
								e.printStackTrace();
							}
							break;
						default:
							System.out.println("function type error ...");
							System.exit(0);
						}
					}
				});
			}
		}
		pool.shutdown();
	}

	@Test
	public void runTest() {
		
		Constants.intializeProject();
		int numOfThreads = 4;
		run(numOfThreads);
	}

}
