package edu.albany.cs.scoreFuncs;

import org.ujmp.core.DenseMatrix;
import org.ujmp.core.Matrix;
import org.ujmp.core.SparseMatrix;

import java.math.BigDecimal;
import java.util.ArrayList;

public class LagrangianFunction implements Function {

	private final Matrix yk;
	private final Matrix zk;
	private final ArrayList<Pair[]> FPairs;
	private final ArrayList<Pair[]> FtPairs;
	private final double rho;
	private final Function func;

	private final FuncType funcID;

	public LagrangianFunction(Function func, Matrix yk, Matrix zk, ArrayList<Pair[]> F, ArrayList<Pair[]> Ft,
			double rho) {
		this.yk = yk;
		this.zk = zk;
		this.FPairs = F;
		this.FtPairs = Ft;
		this.rho = rho;
		this.func = func;
		funcID = FuncType.LagrangianFunc;
	}

	@Override
	public double[] getGradient(double[] x) {

		double[] projectedX = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			projectedX[i] = x[i];
			if (x[i] < 0.0D) {
				projectedX[i] = 0.0D;
			}
			if (x[i] > 1.0D) {
				projectedX[i] = 1.0D;
			}
		}
		double[] gradientFx = func.getGradient(projectedX);
		/**
		 * in lagrangian function, we need to minimize the function, which means
		 * we need to get the negative gradient.
		 */
		Matrix xMat = SparseMatrix.Factory.zeros(x.length, 1);
		Matrix item1 = SparseMatrix.Factory.zeros(x.length, 1);
		for (int i = 0; i < gradientFx.length; i++) {
			gradientFx[i] = -gradientFx[i];
			if (x[i] != 0.0D) {
				xMat.setAsDouble(x[i], i, 0);
			}
			if (gradientFx[i] != 0.0D) {
				item1.setAsDouble(gradientFx[i], i, 0);
			}
		}
		Matrix item2 = multiPlyFtx(FtPairs, yk);
		Matrix Fx = multiPlyFx(FPairs, xMat);
		Matrix FxMinusZk = Fx.minus(zk);
		Matrix FtFxmatZk = multiPlyFtx(FtPairs, FxMinusZk);
		Matrix item3 = FtFxmatZk.times(rho);
		Matrix result = item1.plus(item2).plus(item3);
		double[] resultDouble = new double[x.length];
		int count = 0;
		for (long[] coordinates : result.allCoordinates()) {
			resultDouble[count++] = result.getAsDouble(coordinates);
		}
		return resultDouble;
	}

	private Matrix multiPlyFx(ArrayList<Pair[]> FPairs, Matrix x) {
		Matrix result = SparseMatrix.Factory.zeros(FPairs.size(), 1);
		for (int i = 0; i < FPairs.size(); i++) {
			Pair[] pairs = FPairs.get(i);
			int index0 = pairs[0].i;
			double value0 = pairs[0].value;
			int index1 = pairs[1].i;
			double value1 = pairs[1].value;
			result.setAsDouble(value0 * x.getAsDouble(index0, 0) + value1 * x.getAsDouble(index1, 0), i, 0);
		}
		return result;
	}

	private Matrix multiPlyFtx(ArrayList<Pair[]> FtPairs, Matrix x) {
		Matrix result = SparseMatrix.Factory.zeros(FtPairs.size(), 1);
		for (int i = 0; i < FtPairs.size(); i++) {
			double sum = 0.0D;
			for (Pair pair : FtPairs.get(i)) {
				double value = pair.value;
				int index = pair.i;
				sum += value * x.getAsDouble(index, 0);
			}
			result.setAsDouble(sum, i, 0);
		}
		return result;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see scoreFuncs.Function#getFuncValue(double[]) check the following link
	 * :
	 */
	@Override
	public double getFuncValue(double[] x) {
		// for item1
		double[] dx = new double[x.length];
		for (int i = 0; i < dx.length; i++) {
			dx[i] = x[i];
			if (dx[i] < 0.0D) {
				dx[i] = 0.0D;
			}
			if (dx[i] > 1.0D) {
				dx[i] = 1.0D;
			}
		}
		double item1 = -func.getFuncValue(dx);

		// for item2
		Matrix xMat = SparseMatrix.Factory.zeros(x.length, 1);
		for (int i = 0; i < x.length; i++) {
			xMat.setAsDouble(x[i], i, 0);
		}
		Matrix FPairsxMat = multiPlyFx(FPairs, xMat);
		Matrix FPairsxMatMzk = FPairsxMat.minus(zk);
		double sumItem2 = 0.0D;
		for (int i = 0; i < yk.getRowCount(); i++) {
			sumItem2 += yk.getAsDouble(i, 0) * FPairsxMatMzk.getAsDouble(i, 0);
		}
		// for item3
		double doubleItem3 = 0.0D;
		Matrix item3 = FPairsxMatMzk;
		if (item3.getColumnCount() == 1) {
			for (long[] coordinates : item3.allCoordinates()) {
				doubleItem3 += Math.pow(item3.getAsDouble(coordinates), 2);
			}
		}
		doubleItem3 = (rho / 2.0D) * doubleItem3;
		double lagFuncValue = item1 + sumItem2 + doubleItem3;
		return lagFuncValue;
	}

	@Override
	public double[] getArgMaxFx(ArrayList<Integer> S) {
		return null;
	}

	@Override
	public FuncType getFuncID() {
		return funcID;
	}

	public static void main(String args[]) {
		double[] x = new double[] { 1.0, 2.0, 3.0 };
		Matrix xx = DenseMatrix.Factory.importFromArray(x).transpose();
		System.out.println(xx.toString());
		System.out.println(xx.getRowCount());
		System.out.println(xx.getColumnCount());
	}

	@Override
	public double[] getGradient(int[] S) {
		return null;
	}

	@Override
	public double getFuncValue(int[] S) {
		return 0;
	}

	@Override
	public BigDecimal[] getGradientBigDecimal(BigDecimal[] x) {
		double[] dx = new double[x.length];
		for (int i = 0; i < dx.length; i++) {
			dx[i] = x[i].doubleValue();
			if (dx[i] < 0.0D) {
				dx[i] = 0.0D;
			}
			if (dx[i] > 1.0D) {
				dx[i] = 1.0D;
			}
		}
		double[] gradientFx = func.getGradient(dx);
		/**
		 * in lagrangian function, we need to minimize the function, which means
		 * we need to get the negative gradient.
		 */
		Matrix xMat = SparseMatrix.Factory.zeros(x.length, 1);
		Matrix item1 = SparseMatrix.Factory.zeros(x.length, 1);
		for (int i = 0; i < gradientFx.length; i++) {
			gradientFx[i] = -gradientFx[i];
			if (x[i].doubleValue() != 0.0D) {
				xMat.setAsBigDecimal(x[i], i, 0);
			}
			if (gradientFx[i] != 0.0D) {
				item1.setAsDouble(gradientFx[i], i, 0);
			}
		}
		Matrix item2 = multiPlyFtx(FtPairs, yk);
		Matrix Fx = multiPlyFx(FPairs, xMat);
		Matrix FxMinusZk = Fx.minus(zk);
		Matrix FtFxmatZk = multiPlyFtx(FtPairs, FxMinusZk);
		Matrix item3 = FtFxmatZk.times(rho);
		Matrix result = item1.plus(item2).plus(item3);
		BigDecimal[] resultDouble = new BigDecimal[x.length];
		int count = 0;
		for (long[] coordinates : result.allCoordinates()) {
			resultDouble[count++] = result.getAsBigDecimal(coordinates);
		}
		return resultDouble;
	}

	public class Pair {
	    public int i;
	    public double value;
	    
	    public Pair(int i, double value) {
	        this.i = i;
	        this.value = value;
	    }
	}
}
