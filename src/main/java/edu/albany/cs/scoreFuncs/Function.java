package edu.albany.cs.scoreFuncs;

import java.math.BigDecimal;
import java.util.ArrayList;

/**
 * @author Baojian bzhou6@albany.edu
 */
public interface Function {

	/**
	 * @param x
	 *            get gradient of this function at point x
	 * @return the gradient vector of this function
	 */
	public double[] getGradient(double[] x);

	/**
	 * @param x
	 *            get gradient of this function at point x
	 * @return the gradient vector of this function
	 */
	public BigDecimal[] getGradientBigDecimal(BigDecimal[] x);

	/**
	 * @param S
	 *            get gradient of this function at point S S can be converted to
	 *            x : x[i] = 1.0D if i \in S; x[i] = 0.0D if i \notin S;
	 * @return the gradient vector of this function
	 */
	public double[] getGradient(int[] S);

	/**
	 * @param x
	 *            get value of this function at point x
	 * @return the value of this function
	 */
	public double getFuncValue(double[] x);

	/**
	 * @param S
	 *            get value of this function at set x
	 * @return the value of this function
	 */
	public double getFuncValue(int[] S);

	/**
	 * @param S
	 *            get vector of argmin_x f(x) at point S S can be converted to x
	 *            : x[i] = 1.0D if i \in S; x[i] = 0.0D if i \notin S;
	 * @return the minimizer vector of this function
	 */
	public double[] getArgMaxFx(ArrayList<Integer> S);

	/**
	 * @return the function ID
	 */
	public FuncType getFuncID();
}
