package edu.albany.cs.scoreFuncs;

import edu.albany.cs.base.ArrayIndexComparator;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class KulldorffStat implements Function {

	private final double[] b;
	private final double[] c;
	private final double BAll;
	private final double CAll;
	private final FuncType funcID;
	private final int n;

	private int verboseLevel = 0;

	public KulldorffStat(double[] b, double[] c) {

		this.funcID = FuncType.Kulldorff;
		if (!checkInput(b, c)) {
			System.out.println(funcID + " input parameter is invalid.");
			System.exit(0);
		}
		this.b = b;
		this.c = c;
		this.n = b.length;
		BAll = StatUtils.sum(b);
		CAll = StatUtils.sum(c);
	}

	private boolean checkInput(double[] b, double[] c) {

		if (b == null || c == null) {
			return false;
		} else {
			for (int i = 0; i < n; i++) {
				if (b[i] <= 0.0D) {
					return false;
				}
				if (c[i] < 0.0D) {
					return false;
				}
			}
			return true;
		}
	}

	/**
	 * calculate the gradient of function f denoted in example 2
	 *
	 * @param x
	 *            vector of parameters
	 * @param c
	 *            vector constant values
	 * @param b
	 *            vector constant values
	 * @return the gradient of function f denoted in example 2
	 */
	public double[] getGradient(double[] x) {

		checkIndictorVect(x);
		double[] gradient = new double[n];
		double B = new ArrayRealVector(x).dotProduct(new ArrayRealVector(b));
		double C = new ArrayRealVector(x).dotProduct(new ArrayRealVector(c));
		if (verboseLevel > 0) {
			System.out.println("CAll: " + CAll + " ; BAll: " + BAll + " ; sum of x: " + StatUtils.sum(x));
			System.out.println("B: " + B + "C: " + C);
			System.out.println("(BAll - B): " + (BAll - B) + " ;(CAll-C): " + (CAll - C));
		}
		if ((C / B) > (CAll / BAll)) {
			for (int i = 0; i < n; i++) {
				double item1 = c[i] * (Math.log(C / B) + Math.log((BAll - B) / (CAll - C)));
				double item2 = b[i] * ((CAll - C) / (BAll - B) - C / B);
				gradient[i] = item1 + item2;
				
				if (!Double.isFinite(gradient[i])) {
					System.out.println(gradient[i]);
					System.out.println("gradient error in kulldorff stat");
					System.exit(0);
				}
			}
		} else {
			Arrays.fill(gradient, 0.0D);
			if (verboseLevel > 0) {
				System.out.println("gradient is 0 ...");
			}
		}
		return gradient;
	}

	private void checkIndictorVect(double[] x) {

		if (x == null || x.length != n) {
			System.out.println("Kulldorff gradient error : Invalid parameters ...");
			System.exit(0);
		}
		/** make sure x is an indicator vector */
		for (int i = 0; i < n; i++) {
			if (x[i] < 0.0D || x[i] > 1.0D) {
				System.out.println("x[i] should be in [0,1], but it is " + x[i]);
				System.exit(0);
			}
		}
	}

	/**
	 * get function value of Kulldorff
	 * 
	 */
	public double getFuncValue(double[] x) {

		checkIndictorVect(x);
		double B = new ArrayRealVector(x).dotProduct(new ArrayRealVector(b));
		double C = new ArrayRealVector(x).dotProduct(new ArrayRealVector(c));
		if (verboseLevel > 0) {
			System.out.println("BAll is : " + BAll + " ; CAll is : " + CAll + " ; X is : " + StatUtils.sum(x));
			System.out.println("B is : " + B + " ; C is : " + C);
			System.out.println("BAll - B is : " + (BAll - B) + " ; CAll - C is : " + (CAll - C));
		}
		double item1 = calXlogXA(C, B);
		double item2 = calXlogXA(CAll - C, BAll - B);
		double item3 = calXlogXA(CAll, BAll);
		;
		double f = 0.0D;
		if ((C / B) > (CAll / BAll)) {
			f = item1 + item2 - item3;
			if (verboseLevel > 0) {
				System.out.println("B: " + B);
				System.out.println("C: " + C);
				System.out.println("CAll: " + CAll);
				System.out.println("BAll: " + BAll);
				System.out.println("C / B: " + (C / B));
				System.out.println("CAll / BAll: " + (CAll / BAll));
			}
			if (!Double.isFinite(f)) {
				System.out.println("B: " + B);
				System.out.println("C: " + C);
				System.out.println("CAll: " + CAll);
				System.out.println("BAll: " + BAll);
				System.out.println("Error : Kulldorff value f is " + f);
				System.exit(0);
			}
		} else {
			if (verboseLevel > 0) {
				System.out.println("funcVal something is bad ...");
				System.out.println("C / B: " + C / B);
				System.out.println("CAll / BAll: " + CAll / BAll);
			}
		}
		return f;
	}
	
	public double getRatio(double[] x){
		double B = new ArrayRealVector(x).dotProduct(new ArrayRealVector(b));
		double C = new ArrayRealVector(x).dotProduct(new ArrayRealVector(c));
		return C / B ;
	}

	/**
	 * calculate function xlog(x/a)
	 */
	private double calXlogXA(double x, double a) {
		if (x <= 0.0D) {
			return 0.0D;
		} else {
			if (a <= 0.0D) {
				System.out.println("function xlog(x/a) is error");
				System.exit(0);
			}
			return x * Math.log(x / a);
		}
	}

	/** TODO maximize or minimize */
	public double[] getArgMaxFx(ArrayList<Integer> S) {

		double[] result = new double[n];
		Double[] vectorRatioCB = new Double[S.size()];
		for (int i = 0; i < S.size(); i++) {
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
