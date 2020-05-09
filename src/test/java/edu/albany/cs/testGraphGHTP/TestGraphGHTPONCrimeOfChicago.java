package edu.albany.cs.testGraphGHTP;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.Constants;
import edu.albany.cs.base.Stat;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graphGHTP.GraphGHTP;
import edu.albany.cs.scoreFuncs.FuncType;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.ScoreFuncFactory;
import org.junit.Test;

public class TestGraphGHTPONCrimeOfChicago {

	private final String resultFileName = Constants.ChicagoCrimeOutputFolder + "graph_GHTP_CrimeOfChicago_Result.txt";

	private int verboseLevel = 0;

	private void run(int numOfThreads, FuncType[] funcs) {

		ExecutorService pool = Executors.newFixedThreadPool(numOfThreads);

		for (final File apdmFile : new File(Constants.ChicagoCrimeDataFolder).listFiles()) {

			final APDMInputFormat apdm = new APDMInputFormat(apdmFile);
			final int graphSize = apdm.data.numNodes;
			final ArrayList<Integer[]> edges = apdm.data.intEdges;
			final ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;

			for (final FuncType funcType : funcs) {

				pool.execute(new Thread() {
					public void run() {
						System.out.println("processing : " + apdmFile.getName() + " funcID: " + funcType);
						long startTime = System.nanoTime();
						Function func = null;
						double[] b = null;
						double[] c = null;
						int[] trueSubGraph = apdm.data.trueSubGraphNodes;
						switch (funcType) {
						case Kulldorff:
							b = new double[apdm.data.numNodes];
							c = apdm.data.counts;
							Arrays.fill(b, StatUtils.max(c));
							func = ScoreFuncFactory.getFunc(funcType, b, c);
							
							double bestFuncValue = -Double.MAX_VALUE;
							double[] bestFuncs = null;
							GraphGHTP bestGraphGHTP = null;
							
							for (int s = 50; s <= 1000; s += 50) {
								int g = 1;
								double B = s - 1 + 0.0D;
								boolean isInitialSingle = false;
								GraphGHTP graphGHTP = new GraphGHTP(graphSize, edges, edgeCosts, c, s, g, B,
										isInitialSingle, trueSubGraph, func, true, 1.0D);
								if (bestFuncValue < graphGHTP.funcValueTail) {
									bestFuncValue = graphGHTP.funcValueTail;
									bestGraphGHTP = graphGHTP;
									bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
								}
								if (verboseLevel == 0) {
									System.out.println("-------------------------------------------------");
									System.out.println("processing sparsity s : " + s);
									System.out.println("current function Value is: " + graphGHTP.funcValueTail);
									System.out.println("number of nodes: " + graphGHTP.resultNodesTail.size());
									System.out.println("current best function Value is: " + bestFuncValue);
									System.out.println("-------------------------------------------------");
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
								fileWriter.write(apdmFile.getName() + "," + funcType + "," + funcVal + "," + runTime + ","
										+ pre + "," + rec + "," + fmeasure + "," + fValues + "," + numOfNodes + "\n");
								fileWriter.close();
							} catch (IOException e) {
								e.printStackTrace();
							}
							break;

						case EMS:
							b = new double[apdm.data.numNodes];
							c = new double[apdm.data.numNodes];
							Arrays.fill(b, 1.0D);
							double max = StatUtils.max(apdm.data.counts);
							for (int i = 0; i < apdm.data.numNodes; i++) {
								c[i] = apdm.data.counts[i] / max;
							}
							func = ScoreFuncFactory.getFunc(funcType, b, c);
							
							bestFuncValue = -Double.MAX_VALUE;
							bestFuncs = null;
							bestGraphGHTP = null;
							
							for (int s = 50; s <= 1000; s += 50) {
								int g = 1;
								double B = s - 1 + 0.0D;
								boolean isInitialSingle = false;
								GraphGHTP graphGHTP = new GraphGHTP(graphSize, edges, edgeCosts, c, s, g, B,
										isInitialSingle, trueSubGraph, func, true, 1.0D);
								if (bestFuncValue < graphGHTP.funcValueTail) {
									bestFuncValue = graphGHTP.funcValueTail;
									bestGraphGHTP = graphGHTP;
									bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
								}
								if (verboseLevel == 0) {
									System.out.println("-------------------------------------------------");
									System.out.println("processing sparsity s : " + s);
									System.out.println("current function Value is: " + graphGHTP.funcValueTail);
									System.out.println("number of nodes: " + graphGHTP.resultNodesTail.size());
									System.out.println("current best function Value is: " + bestFuncValue);
									System.out.println("-------------------------------------------------");
								}
							}
							try {
								double runTime = (System.nanoTime() - startTime) / 1e9;
								double funcVal = bestGraphGHTP.funcValueTail;
								double pre = 0.0D;
								double rec = 0.0D;
								double fmeasure = 0.0D;
								String fValues = "";
								int numOfNodes = bestGraphGHTP.resultNodesTail.size();
								for (double fVal : bestGraphGHTP.fValues) {
									fValues += fVal + " ";
								}
								FileWriter fileWriter = new FileWriter(new File(resultFileName), true);
								fileWriter.write(apdmFile.getName() + "," + funcType + "," + funcVal + "," + runTime + ","
										+ pre + "," + rec + "," + fmeasure + "," + fValues + "," + numOfNodes + "\n");
								fileWriter.close();
							} catch (IOException e) {
								e.printStackTrace();
							}
							break;

						case EBP:
							b = apdm.data.averCounts;
							c = apdm.data.counts;
							double[] nonZeros = new Stat(apdm.data.counts).nonZeros();
							Arrays.fill(b, StatUtils.mean(nonZeros));
							System.out.println(StatUtils.mean(c));
							Arrays.fill(b,StatUtils.mean(c));
							func = ScoreFuncFactory.getFunc(funcType, b, c);
							
							bestFuncValue = -Double.MAX_VALUE;
							bestFuncs = null;
							bestGraphGHTP = null;

							for (int s = 50; s <= 1000; s += 50) {
								int g = 1;
								double B = s - 1 + 0.0D;
								boolean isInitialSingle = false;
								GraphGHTP graphGHTP = new GraphGHTP(graphSize, edges, edgeCosts, c, s, g, B,
										isInitialSingle, trueSubGraph, func, true, 1.0D);
								if (bestFuncValue < graphGHTP.funcValueTail) {
									bestFuncValue = graphGHTP.funcValueTail;
									bestGraphGHTP = graphGHTP;
									bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
								}
								if (verboseLevel == 0) {
									System.out.println("-------------------------------------------------");
									System.out.println("processing sparsity s : " + s);
									System.out.println("current function Value is: " + graphGHTP.funcValueTail);
									System.out.println("number of nodes: " + graphGHTP.resultNodesTail.size());
									System.out.println("current best function Value is: " + bestFuncValue);
									System.out.println("-------------------------------------------------");
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
								fileWriter.write(apdmFile.getName() + "," + funcType + "," + funcVal + "," + runTime + ","
										+ pre + "," + rec + "," + fmeasure + "," + fValues + "," + numOfNodes + "\n");
								fileWriter.close();
							} catch (IOException e) {
								e.printStackTrace();
							}
							System.exit(0);
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

	public static void test() {
		APDMInputFormat apdm = new APDMInputFormat("data/CrimeOfChicago/graph/APDM-2015.txt");
		int count = 0;
		double sumAver = 0;
		double sumCounts = 0;
		int count2 = 0;
		double sumAver2 = 0;
		double sumCounts2 = 0;
		double nonZeros = 0.0D;
		int index = 0;
		for (int i = 0; i < apdm.data.numNodes; i++) {
			if (apdm.data.counts[i] != 0.0D) {
				nonZeros += apdm.data.counts[i];
				index++;
			}
		}
		Arrays.fill(apdm.data.averCounts, nonZeros / (index + 0.0D));
		for (int i = 0; i < apdm.data.numNodes; i++) {
			if (apdm.data.averCounts[i] > apdm.data.counts[i] && apdm.data.averCounts[i] > 1.0D) {
				System.out.println(apdm.data.averCounts[i] + " " + apdm.data.counts[i]);
				count++;
				sumAver += apdm.data.averCounts[i];
				sumCounts += apdm.data.counts[i];
			}
			if (apdm.data.averCounts[i] < apdm.data.counts[i]) {
				System.out.println(apdm.data.averCounts[i] + " " + apdm.data.counts[i]);
				count2++;
				sumAver2 += apdm.data.averCounts[i];
				sumCounts2 += apdm.data.counts[i];
			}
		}
		System.out.println("total count is: " + count);
		System.out.println("sumAver is: " + sumAver);
		System.out.println("sumCounts is: " + sumCounts);

		System.out.println("total count2 is: " + count2);
		System.out.println("sumAver2 is: " + sumAver2);
		System.out.println("sumCounts2 is: " + sumCounts2);
		Utils.stop();
	}

	@Test
	public void runTest() {
		
		Constants.intializeProject();
		int numOfThreads = 1;
		run(numOfThreads, new FuncType[] { FuncType.EBP, FuncType.EMS, FuncType.Kulldorff, });
	}

}
