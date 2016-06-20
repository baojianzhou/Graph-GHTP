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
import edu.albany.cs.graphIHT.GraphIHT;
import edu.albany.cs.scoreFuncs.FuncType;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.ScoreFuncFactory;

public class TestGraphIHTOnCitHep {

	private final int numOfThreads;
	private final FuncType[] funcs;

	private final String citHepPhFolder = "data/CitHepPh/testGraph";
	private final String resultFileName = "./output/CitHepPh/graph_IHT_CitHepPh_Single.txt";

	private int verboseLevel = 0;

	public TestGraphIHTOnCitHep(int numOfThreads) {
		this.numOfThreads = numOfThreads;
		this.funcs = new FuncType[] { FuncType.EMS };
		run();
	}

	private void run() {

		if (!new File("./output/CitHepPh/").isDirectory()) {
			new File("./output/CitHepPh/").mkdir();
		}
		ExecutorService pool = Executors.newFixedThreadPool(numOfThreads);

		for (final File apdmFile : new File(citHepPhFolder).listFiles()) {
			
			final APDMInputFormat apdm = new APDMInputFormat(apdmFile);
			final int graphSize = apdm.data.numNodes;
			final ArrayList<Integer[]> edges = apdm.data.intEdges;
			final ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;

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
							b = apdm.data.averCounts;
							c = apdm.data.counts;
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
							System.out.println("BAll: " + StatUtils.sum(b));
							System.out.println("CAll: " + StatUtils.sum(c));
						}

						double bestFuncValue = -Double.MAX_VALUE;
						double[] bestFuncs = null;
						GraphIHT bestGraphIHT = null;

						for (int s = 50; s <= 1000; s += 50) {
							int g = 1;
							double B = s - 1 + 0.0D;
							GraphIHT graphIHT = new GraphIHT(graphSize, edges, edgeCosts, apdm.data.counts, s, g, B,
									false, trueSubGraph, func);
							System.out.println("fValues: "+graphIHT.fValues.toString());
							if (bestFuncValue < graphIHT.funcValueTail) {
								bestFuncValue = graphIHT.funcValueTail;
								bestGraphIHT = graphIHT;
								bestFuncs = ArrayUtils.add(bestFuncs, bestFuncValue);
							} else {
							}
						}

						try {
							double runTime = (System.nanoTime() - startTime) / 1e9;
							double funcVal = bestGraphIHT.funcValueTail;
							double pre = 0.0D;
							double rec = 0.0D;
							double fmeasure = 0.0D;
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
		for (File file : new File(citHepPhFolder).listFiles()) {
			allFilePaths = ArrayUtils.add(allFilePaths, citHepPhFolder + "/" + file.getName());
		}
		return allFilePaths;
	}

	public static void main(String args[]) {
		if (args == null || args.length == 0) {
			int numOfThreads = 3;
			new TestGraphIHTOnCitHep(numOfThreads);
		} else {
			int numOfThreads = Integer.parseInt(args[0]);
			new TestGraphIHTOnCitHep(numOfThreads);
		}
	}

}
