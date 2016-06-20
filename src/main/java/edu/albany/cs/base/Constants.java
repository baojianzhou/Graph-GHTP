package edu.albany.cs.base;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import edu.albany.cs.scoreFuncs.FuncType;

/**
 * Some constants in this project
 * 
 * @author baojian bzhou6@albany.edu
 *
 */
public final class Constants {

	private Constants() {
	}

	public static final String line = "----------------------------------------------------------";
	public static final String algo1 = "Graph-IHT";
	public static final String algo2 = "Graph-GHTP";

	/** data folder */
	public static final String BWSNDataFolder = "./data/BWSN/testGraphs";
	public static final String CitHepPhDataFolder = "./data/CitHepPh/testGraphs";
	public static final String TrafficDataFolder = "./data/Traffic/testGraphs";
	public static final String ChicagoCrimeDataFolder = "./data/ChicagoCrime/testGraphs";

	/** result folder */
	public static final String output = "./output";
	public static final String BWSNOutputFolder = output + "/BWSN/";
	public static final String CitHepPhOutputFolder = output + "/CitHepPh/";
	public static final String TrafficOutputFolder = output + "/Traffic/";
	public static final String ChicagoCrimeOutputFolder = output + "/ChicagoCrime/";

	/** score functions */
	public static final FuncType[] funcs = new FuncType[] { FuncType.Kulldorff, FuncType.EMS, FuncType.EBP };

	
	
	/**
	 * 1. check the data folder
	 * 2. create result folders if they haven't been created.
	 */
	public static void intializeProject() {
		System.out.println(line);
		System.out.println("initialize project ....");
		if (Files.exists(Paths.get(BWSNDataFolder))) {
			System.out.println("BWSN data location: " + BWSNDataFolder);
			for (String file : new File(BWSNDataFolder).list()) {
				System.out.println("file: " + file);
			}
		} else {
			exit("cannot find BWSN data folder", 0);
		}
		if (Files.exists(Paths.get(CitHepPhDataFolder))) {
			System.out.println("CitHepPh data location: " + CitHepPhDataFolder);
			for (String file : new File(CitHepPhDataFolder).list()) {
				System.out.println("file: " + file);
			}
		} else {
			exit("cannot find CitHepPh data folder", 0);
		}
		if (Files.exists(Paths.get(TrafficDataFolder))) {
			System.out.println("Traffic data location: " + TrafficDataFolder);
			for (String file : new File(TrafficDataFolder).list()) {
				System.out.println("file: " + file);
			}
		} else {
			exit("cannot find CitHepPh data folder", 0);
		}
		if (Files.exists(Paths.get(ChicagoCrimeDataFolder))) {
			System.out.println("ChicagoCrime data location: " + ChicagoCrimeDataFolder);
			for (String file : new File(ChicagoCrimeDataFolder).list()) {
				System.out.println("file: " + file);
			}
		} else {
			exit("cannot find CitHepPh data folder", 0);
		}
		new File(BWSNOutputFolder).mkdirs();
		new File(CitHepPhOutputFolder).mkdirs();
		new File(TrafficOutputFolder).mkdirs();
		new File(ChicagoCrimeOutputFolder).mkdirs();
	}

	public static void exit(String errorInfo, int status) {
		System.out.println("error Info: " + errorInfo);
		System.exit(status);
	}

	public static void main(String args[]) {
		intializeProject();
	}

}
