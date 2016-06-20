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
import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graphGHTP.GraphGHTP;
import edu.albany.cs.scoreFuncs.FuncType;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.ScoreFuncFactory;

public class TestGraphGHTPONBWSN {
	private final int numOfThreads;
	private final FuncType[] funcs;
	private final String waterNetworkDataFolder = "data/BWSN/source_12500";
	private final String resultFileName = "./output/BWSN/graph_GHTP_BWSN_Result.txt";

	private int verboseLevel = 0;

	public TestGraphGHTPONBWSN(int numOfThreads) {
		this.numOfThreads = numOfThreads;
		this.funcs = new FuncType[] { FuncType.Kulldorff, FuncType.EMS, FuncType.EBP };
		run();
	}

	private void run() {

		if (!new File("./output/BWSN/").isDirectory()) {
			new File("./output/BWSN/").mkdir();
		}
		ExecutorService pool = Executors.newFixedThreadPool(numOfThreads);

		for (final File apdmFile : new File(waterNetworkDataFolder).listFiles()) {

			final APDMInputFormat apdm = new APDMInputFormat(apdmFile);
			final int graphSize = apdm.data.numNodes;
			final ArrayList<Integer[]> edges = apdm.data.intEdges;
			final ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;
			if (verboseLevel > 0) {
				if (!apdmFile.getName().equals("APDM-Water-source-12500_time_03_hour_noise_4.txt")) {
					continue;
				}
			}
			final double noiseLevel = Double.parseDouble(apdmFile.getName().split("_")[5].split(".txt")[0]) / 100.0D;
			for (final FuncType funcID : funcs) {

				pool.execute(new Thread() {
					public void run() {
						System.out.println("processing file: " + apdmFile.getName());
						long startTime = System.nanoTime();
						Function func = null;
						double[] b = new double[apdm.data.numNodes];
						double[] c = apdm.data.counts;
						int[] trueSubGraph = apdm.data.trueSubGraphNodes;
						switch (funcID) {
						case Kulldorff:
							Arrays.fill(b, 1.0D);
							func = ScoreFuncFactory.getFunc(funcID, b, c);
							break;
						case EMS:
						case EBP:
							Arrays.fill(b, 0.5D);
							func = ScoreFuncFactory.getFunc(funcID, b, c);
							break;
						default:
							System.out.println("function type error ...");
							System.exit(0);
						}

						if (verboseLevel > 0) {
							System.out.println("---------");
							System.out.println("processing file: " + apdmFile.getName() + " Func: " + funcID);
						}
						if (verboseLevel > 0) {
							System.out.println("noiseLevel: " + noiseLevel);
							System.out.println("BAll: " + StatUtils.sum(b));
							System.out.println("CAll: " + StatUtils.sum(c));
						}

						int[] S = null;
						for (int i = 50; i <= 1000; i += 50) {
							S = ArrayUtils.add(S, i);
						}
						double bestFuncValue = -Double.MAX_VALUE;
						double[] bestFuncs = null;
						GraphGHTP bestGraphGHTP = null;
						for (int s : S) {
							int g = 1;
							double B = s - 1 + 0.0D;
							boolean isInitialSingle = false;
							GraphGHTP graphGHTP = new GraphGHTP(graphSize, edges, edgeCosts, apdm.data.counts, s, g, B,
									isInitialSingle, trueSubGraph, func, true, 1.0D);
							if (bestFuncValue < graphGHTP.funcValueTail) {
								bestFuncValue = graphGHTP.funcValueTail;
								bestGraphGHTP = graphGHTP;
								bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
							}
							System.out.println("s: " + s);
							if (verboseLevel > 0) {
								System.out.println("-------------------------------------------------");
								PreRec preRec = new PreRec(graphGHTP.resultNodesTail,trueSubGraph);
								System.out.println(preRec.toString());
								System.out.println("current function Value is: " + graphGHTP.funcValueTail);
								System.out.println("current best function Value is: " + bestFuncValue);
								System.out.println("running Time: " + graphGHTP.runTime);
								System.out.println("-------------------------------------------------");
								Utils.stop();
							}
						}
						if (verboseLevel == 0) {
							System.out.println("-------------------------------------------------");
							System.out.println("fValues: " + bestGraphGHTP.fValues.toString());
							System.out.println("bestFuncs: " + bestFuncs);
							System.out.println("running Time: " + bestGraphGHTP.runTime);
							System.out.println("length of result: " + bestGraphGHTP.resultNodesTail.size());
							System.out.println("length of trueSubGraph: " + apdm.data.trueSubGraphNodes.length);
							PreRec preRec = new PreRec(bestGraphGHTP.resultNodesTail, trueSubGraph);
							System.out.println(preRec.toString());
							System.out.println("-------------------------------------------------");

							double bb = 0.0D;
							double cc = 0.0D;
							for (int i : apdm.data.trueSubGraphNodes) {
								bb += b[i];
							}
							for (int i : apdm.data.trueSubGraphNodes) {
								cc += c[i];
							}
							System.out.println("B: " + bb + " C: " + cc + " C / B: " + (cc / bb));

							bb = 0.0D;
							cc = 0.0D;
							for (int i : bestGraphGHTP.resultNodesTail) {
								bb += b[i];
							}
							for (int i : bestGraphGHTP.resultNodesTail) {
								cc += c[i];
							}
							System.out.println("B: " + bb + " C: " + cc + " C / B: " + (cc / bb));
						}

						try {
							double runTime = (System.nanoTime() - startTime) / 1e9;
							PreRec preRec = new PreRec(bestGraphGHTP.resultNodesTail, trueSubGraph);
							double funcVal = bestGraphGHTP.funcValueTail;
							double pre = preRec.pre;
							double rec = preRec.rec;
							double fmeasure = preRec.fmeasure;
							String fValues = "";
							for (double fVal : bestGraphGHTP.fValues) {
								fValues += fVal + " ";
							}
							FileWriter fileWriter = new FileWriter(new File(resultFileName), true);
							fileWriter.write(apdmFile.getName() + "," + funcID + "," + funcVal + "," + runTime + ","
									+ pre + "," + rec + "," + fmeasure + "," + fValues + "\n");
							fileWriter.close();
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				});
			}
		}
		pool.shutdown();
	}

	public static void main(String args[]) {
		if (args == null || args.length == 0) {
			int numOfThreads = 5;
			new TestGraphGHTPONBWSN(numOfThreads);
		} else {
			int numOfThreads = Integer.parseInt(args[0]);
			new TestGraphGHTPONBWSN(numOfThreads);
		}
	}

}
