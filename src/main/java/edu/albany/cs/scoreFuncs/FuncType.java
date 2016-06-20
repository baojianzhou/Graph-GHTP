package edu.albany.cs.scoreFuncs;

/**
 * Score functions that we have now.
 *
 * @author baojian bzhou6@albany.edu
 */
public enum FuncType {

	/** Kulldorff Scan Statistic */
	Kulldorff,
	/** Expectation Based Poisson */
	EBP,
	/** Elevated mean statistic */
	EMS,
	/** Lagrangian Function */
	LagrangianFunc,
	/** Unknown Type of Function */
	Unknown;

	public static FuncType defaultFuncType() {
		return FuncType.Unknown;
	}

}
