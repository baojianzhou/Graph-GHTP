package edu.albany.cs.testScoreFuncs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import edu.albany.cs.base.ArrayIndexComparator;
import edu.albany.cs.scoreFuncs.EBPStat;

public class TestEBPStat {
	

	@Test
	public void testGetArgMaxFx() {
		double[] b = new double[] { 1.0D };
		double[] c = new double[] { 1.0D };
		EBPStat ebp = new EBPStat(b, c);
		ArrayList<Integer> S = new ArrayList<Integer>();
		S.add(1);
		S.add(2);
		S.add(3);
		S.add(4);
		double[] result = new double[4];
		Double[] vectorRatioCB = new Double[S.size()];
		vectorRatioCB = new Double[] { 0.2, 0.4, 0.1, 0.3 };
		ArrayIndexComparator arrayIndexComparator = new ArrayIndexComparator(vectorRatioCB);
		Integer[] indexes = arrayIndexComparator.indexes;
		Arrays.sort(indexes, arrayIndexComparator);
		ArrayList<Integer> sortedS = new ArrayList<Integer>(); // v_1,v_2,...,v_m
		for (int index : indexes) {
			sortedS.add(S.get(index));
		}
		System.out.println("indexes: " + Arrays.toString(indexes));
		System.out.println("vectorRatioCB: " + Arrays.toString(vectorRatioCB));
		System.out.println("soredS: " + sortedS.toString());
		double maxF = -Double.MAX_VALUE;
		double[] argMaxX = null;
		for (int k = 1; k <= sortedS.size(); k++) {
			List<Integer> Rk = sortedS.subList(0, k);
			double[] x = new double[ebp.getSize()];
			for (int i = 0; i < ebp.getSize(); i++) {
				x[i] = 0.0D;
			}
			for (int index : Rk) {
				x[index] = 1.0D;
			}
			double fk = ebp.getFuncValue(x);
			if (fk > maxF) {
				maxF = fk;
				argMaxX = x;
			}
		}
		result = argMaxX;
		System.out.println("result: "+Arrays.toString(result));
	}
	
}
