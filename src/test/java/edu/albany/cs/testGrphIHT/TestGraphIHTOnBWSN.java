package edu.albany.cs.testGrphIHT;

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
import edu.albany.cs.graphIHT.GraphIHT;
import edu.albany.cs.scoreFuncs.FuncType;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.ScoreFuncFactory;

public class TestGraphIHTOnBWSN {
	
	//the base of kulldorff is the maximal value of citations

	private final int numOfThreads;
	private final FuncType[] funcs;

	private final String BWSNDataFolder = "data/BWSN/source_12500";
	private final String resultFileName = "./output/BWSN/graph_IHT_BWSN_Result_Updated.txt";

	private int verboseLevel = 0;

	public TestGraphIHTOnBWSN(int numOfThreads) {
		this.numOfThreads = numOfThreads;
		this.funcs = new FuncType[] { FuncType.Kulldorff};
		run();
	}

	private void run() {

		if (!new File("./output/BWSN/").isDirectory()) {
			new File("./output/BWSN/").mkdir();
		}
		ExecutorService pool = Executors.newFixedThreadPool(numOfThreads);

		for (final File apdmFile : new File(BWSNDataFolder).listFiles()) {

			final APDMInputFormat apdm = new APDMInputFormat(apdmFile);
			final int graphSize = apdm.data.numNodes;
			final ArrayList<Integer[]> edges = apdm.data.intEdges;
			final ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;
			final double noiseLevel = Double.parseDouble(apdmFile.getName().split("_")[5].split(".txt")[0]) / 100.0D;

			for (final FuncType funcID : funcs) {

				pool.execute(new Thread() {
					public void run() {
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

						if (verboseLevel == 0) {
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
						GraphIHT bestGraphIHT = null;

						for (int s : S) {
							int g = 1;
							double B = s - 1 + 0.0D;
							GraphIHT graphIHT = new GraphIHT(graphSize, edges, edgeCosts, apdm.data.counts, s, g, B,
									true, trueSubGraph, func);
							if (bestFuncValue < graphIHT.funcValueTail) {
								bestFuncValue = graphIHT.funcValueTail;
								bestGraphIHT = graphIHT;
								bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
								if (verboseLevel == 0) {
									Utils.calPreAndRec(Utils.toArray(graphIHT.resultNodesTail), trueSubGraph);
									System.out.println("function Value is: " + bestFuncValue);
								}
							} else {
								if (verboseLevel == 0) {
									Utils.calPreAndRec(Utils.toArray(graphIHT.resultNodesTail), trueSubGraph);
									System.out.println("function Value is: " + bestFuncValue);
								}
							}
						}

						if (verboseLevel == 0) {
							System.out.println("fValues: " + bestGraphIHT.fValues.toString());
							System.out.println("bestFuncs: " + bestFuncs);
							System.out.println("running Time: " + bestGraphIHT.runTime);
							PreRec preRec = new PreRec(bestGraphIHT.resultNodesTail, trueSubGraph);
							System.out.println(preRec.toString());
						}

						try {
							double runTime = (System.nanoTime() - startTime) / 1e9;
							PreRec preRec = new PreRec(bestGraphIHT.resultNodesTail, trueSubGraph);
							double funcVal = bestGraphIHT.funcValueTail;
							double pre = preRec.pre;
							double rec = preRec.rec;
							double fmeasure = preRec.fmeasure;
							String fValues = "";
							for (double fVal : bestGraphIHT.fValues) {
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

	public String[] getBWSNFilePathList() {
		String[] allFilePaths = null;
		for (File file : new File(BWSNDataFolder).listFiles()) {
			allFilePaths = ArrayUtils.add(allFilePaths, BWSNDataFolder + "/" + file.getName());
		}
		return allFilePaths;
	}

	public static void main(String args[]) {
		if (args == null || args.length == 0) {
			int numOfThreads = 5;
			new TestGraphIHTOnBWSN(numOfThreads);
		} else {
			int numOfThreads = Integer.parseInt(args[0]);
			new TestGraphIHTOnBWSN(numOfThreads);
		}
	}
}
