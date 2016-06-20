package edu.albany.cs.scoreFuncs;

import edu.albany.cs.base.ArrayIndexComparator;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class EMSStat implements Function {

	private final double[] b;
	/** attribute 1 : each node i has a feature c_i */
	private final double[] c;
	private final FuncType funcID;
	private final int n;

	public EMSStat(double[] b, double[] c) {
		funcID = FuncType.EMS;
		if (!checkInput(b, c)) {
			System.out.println(funcID + " input parameter is invalid.");
		}
		this.b = b;
		this.c = c;
		this.n = b.length;
	}

	private boolean checkInput(double[] b, double c[]) {

		if (b == null || c == null) {
			return false;
		} else {
			for (int i = 0; i < b.length; i++) {
				if (b[i] <= 0.0D) {
					return false;
				}
			}
			return true;
		}
	}

	@Override
	public double[] getGradient(double[] x) {

		if (x == null || c == null || x.length != c.length) {
			new IllegalArgumentException("Error : Invalid parameters ...");
			System.exit(0);
		}
		double[] gradient = new double[n];
		double sigmaX = StatUtils.sum(x);
		if (sigmaX == 0.0D) {
			System.out.println(funcID + " Error : the denominator should not be zero.");
			System.exit(0);
		}
		double sigmaCX = new ArrayRealVector(x).dotProduct(new ArrayRealVector(c));
		for (int i = 0; i < gradient.length; i++) {
			// gradient[i] = c[i] * (Math.sqrt(sigmaX) / sigmaX) - (0.5D) *
			// (sigmaCX / Math.pow(sigmaX, 1.5D));
			gradient[i] = c[i] * (Math.sqrt(sigmaX) / sigmaX) - (0.5D) * (sigmaCX / Math.pow(sigmaX, 1.5D));
		}
		return gradient;
	}

	@Override
	public double getFuncValue(double[] x) {

		double funcValue = 0.0D;
		if (x == null || c == null || x.length != c.length) {
			new IllegalArgumentException("Error : Invalid parameters ...");
			System.exit(0);
		}
		double sigmaX = StatUtils.sum(x);
		double sigmaCX = new ArrayRealVector(x).dotProduct(new ArrayRealVector(c));
		if (sigmaX <= 0.0D) {
			System.out.println("funcValue error ...");
			System.exit(0);
		} else {
			funcValue = sigmaCX / Math.sqrt(sigmaX);
		}
		if (!Double.isFinite(funcValue)) {
			System.out.println(funcID + " Error : elevated mean scan stat is not a real value, f is " + funcValue);
			System.exit(0);
		}

		return funcValue;
	}

	@Override
	public double[] getArgMaxFx(ArrayList<Integer> S) {
		double[] result = new double[this.b.length];
		Double[] vectorRatioCB = new Double[S.size()];
		for (int i = 0; i < vectorRatioCB.length; i++) {
			vectorRatioCB[i] = c[S.get(i)] / b[S.get(i)];
		}
		ArrayIndexComparator arrayIndexComparator = new ArrayIndexComparator(vectorRatioCB);
		Integer[] indexes = arrayIndexComparator.indexes;
		Arrays.sort(indexes, arrayIndexComparator);
		ArrayList<Integer> sortedS = new ArrayList<Integer>(); // v_1,v_2,...,v_m
		for (int index : indexes) {
			sortedS.add(S.get(index));
		}
		double maxF = -Double.MAX_VALUE;
		double[] argMaxX = null;
		for (int k = 1; k <= sortedS.size(); k++) {
			List<Integer> Rk = sortedS.subList(0, k);
			double[] x = new double[n];
			for (int i = 0; i < n; i++) {
				x[i] = 0.0D;
			}
			for (int index : Rk) {
				x[index] = 1.0D;
			}
			double fk = getFuncValue(x);
			if (fk > maxF) {
				maxF = fk;
				argMaxX = x;
			}
		}
		result = argMaxX;
		return result;
	}

	@Override
	public FuncType getFuncID() {
		return funcID;
	}

	@Override
	public double[] getGradient(int[] S) {
		double[] x = new double[n];
		Arrays.fill(x, 0.0D);
		for (int i : S) {
			x[i] = 1.0D;
		}
		return getGradient(x);
	}

	@Override
	public double getFuncValue(int[] S) {
		double[] x = new double[n];
		Arrays.fill(x, 0.0D);
		for (int i : S) {
			x[i] = 1.0D;
		}
		return getFuncValue(x);
	}

	@Override
	public BigDecimal[] getGradientBigDecimal(BigDecimal[] x) {
		double[] xD = new double[n];
		for (int i = 0; i < n; i++) {
			xD[i] = x[i].doubleValue();
		}
		double[] gradient = getGradient(xD);
		BigDecimal[] grad = new BigDecimal[n];
		for (int i = 0; i < n; i++) {
			grad[i] = new BigDecimal(gradient[i]);
		}
		return grad;
	}
}
