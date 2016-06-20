package edu.albany.cs.scoreFuncs;

import edu.albany.cs.base.ArrayIndexComparator;
import edu.albany.cs.base.Utils;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * 
 * please see the following paper :
 * http://www.cs.cmu.edu/~./neill/papers/jrssb2012.pdf Fast subset scan for
 * spatial pattern detection The function F(S) is Poisson statistic in table 1
 * 
 * @author baojian bzhou6@albany.edu
 *
 */
public class EBPStat implements Function {

	private final double[] b;
	private final double[] c;
	private final int n;
	private int verboseLevel = 0;

	private final FuncType funcID;

	/**
	 * @param b
	 *            base
	 * @param c
	 *            data (e.g. counts)
	 */
	public EBPStat(double[] b, double[] c) {

		this.funcID = FuncType.EBP;
		if (!checkInput(b, c)) {
			Utils.error(funcID + " input parameter is invalid.", 0);
		}
		this.b = b;
		this.c = c;
		this.n = b.length;

	}

	/**
	 * check EBP statistic inputs is valid
	 * 
	 * @return
	 */
	private boolean checkInput(double[] b, double[] c) {

		if (verboseLevel > 0) {
			System.out.println("BAll: " + StatUtils.sum(b));
			System.out.println("CAll: " + StatUtils.sum(c));
		}
		if (b == null || c == null || b.length == 0 || c.length == 0) {
			return false;
		} else {
			if (StatUtils.sum(b) <= 0.0D || StatUtils.sum(c) <= 0.0D) {
				return false;
			}
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
	 * 
	 */
	@Override
	public double[] getGradient(double[] x) {

		checkIndictorVect(x);
		if (StatUtils.sum(x) <= 0.0D) {
			x[new Random().nextInt(n)] = 1.0D;
		}
		double[] gradient = new double[n];
		double B = new ArrayRealVector(x).dotProduct(new ArrayRealVector(b));
		double C = new ArrayRealVector(x).dotProduct(new ArrayRealVector(c));
		if (verboseLevel > 0) {
			System.out.println("B is : " + B + " ; C is : " + C);
		}
		if (B <= 0.0D || C <= 0.0D) {
			Utils.error("EBP Error : B or C is non-positive value ...\n" + "B is: " + B + " C is: " + C, 0);
		}
		if (C > B) {
			for (int i = 0; i < n; i++) {
				gradient[i] = c[i] * Math.log(C / B) + b[i] * (1 - C / B);
				if (!Double.isFinite(gradient[i])) {
					Utils.error("EBP gradient error ...", 0);
				}
			}
		} else {
			Arrays.fill(gradient, 0.0D);
		}
		return gradient;
	}

	/**
	 * please see the following paper :
	 * http://www.cs.cmu.edu/~./neill/papers/jrssb2012.pdf Fast subset scan for
	 * spatial pattern detection The function F(S) is Poisson statistic in table
	 * 1 vector x is indicator vector x(S) = 1.0 when i \in S, x(\bar{S}) = 0.0
	 */
	@Override
	public double getFuncValue(double[] x) {

		checkIndictorVect(x);
		double B = new ArrayRealVector(x).dotProduct(new ArrayRealVector(b));
		double C = new ArrayRealVector(x).dotProduct(new ArrayRealVector(c));
		double f = 0.0D;
		/** in this situation, the function value is 0.0 */
		if (C > B) {
			f = C * Math.log(C / B) + B - C;
			if (!Double.isFinite(f)) {
				Utils.error("EBP stat getFuncValue error. it is not finite.\n" + "B: " + B + " C: " + C, 0);
			}
		}
		return f;
	}

	private void checkIndictorVect(double[] x) {
		if (x == null || x.length != n) {
			Utils.error("Kulldorff gradient error : Invalid parameters ...", 0);
		}
		/** make sure x is an indicator vector */
		for (int i = 0; i < n; i++) {
			if (x[i] < 0.0D || x[i] > 1.0D) {
				Utils.error("x[i] should be in [0,1], but it is " + x[i], 0);
			}
		}
	}

	/**
	 * To do maximization
	 */
	@Override
	public double[] getArgMaxFx(ArrayList<Integer> S) {

		double[] result = new double[n];
		Double[] vectorRatioCB = new Double[S.size()];
		for (int i = 0; i < S.size(); i++) {
			vectorRatioCB[i] = c[S.get(i)] / b[S.get(i)];
		}
		ArrayIndexComparator arrayIndexComparator = new ArrayIndexComparator(vectorRatioCB);
		Integer[] indexes = arrayIndexComparator.indexes;
		Arrays.sort(indexes, arrayIndexComparator);
		/** v_1,v_2,...,v_m from large to small */
		ArrayList<Integer> sortedS = new ArrayList<Integer>();
		/** indexes from large to small */
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

	public int getSize() {
		return n;
	}

}
