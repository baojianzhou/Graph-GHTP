package edu.albany.cs.scoreFuncs;

/**
 * score function factory for EBP EMS and Kulldorff statistics
 * 
 * @author baojian
 *
 */
public class ScoreFuncFactory {

	public static Function getFunc(FuncType funcID, double[] b, double[] c) {

		if (funcID == null) {
			return null;
		} else if (funcID.equals(FuncType.Kulldorff)) {
			return new KulldorffStat(b, c);
		} else if (funcID.equals(FuncType.EMS)) {
			return new EMSStat(b, c);
		} else if (funcID.equals(FuncType.EBP)) {
			return new EBPStat(b, c);
		} else {
			System.out.println("Unknown Type ...");
			return null;
		}
	}
}
