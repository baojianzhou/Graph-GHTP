package edu.albany.cs.fastPCST;

import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.headApprox.HeadApprox.F;

import java.util.ArrayList;

/**
 * FastPCSF this algorithm is based on FastPCST algorithm
 *
 * @author baojian bzhou6@albany.edu
 */
public class FastPCSF {

    public F f;
    public ArrayList<Integer[]> subTreeEdges = new ArrayList<Integer[]>();
    public ArrayList<Double> subTreeCosts = new ArrayList<Double>();


    public FastPCSF(ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, ArrayList<Double> pi, int g) {

        FastPCST pcstFast = new FastPCST(edges, pi, edgeCosts, FastPCST.kNoRoot, g, FastPCST.PruningMethod.kStrongPruning, -1);
        ArrayList<Integer> resultNodes = null;
        ArrayList<Integer> resultEdges = null;
        if (!pcstFast.run(resultEdges, resultNodes)) {
            System.out.println("Error : Algorithm returned false. There must be an error. \n");
            System.exit(0);
        } else {
            resultEdges = pcstFast.resultEdges; //index
            resultNodes = pcstFast.resultNodes;
        }
        ArrayList<Integer[]> subTreeEdges = new ArrayList<Integer[]>();
        ArrayList<Double> subTreeCosts = new ArrayList<Double>();
        for (int i : pcstFast.resultEdges) {
            subTreeEdges.add(new Integer[]{edges.get(i)[0], edges.get(i)[1]});
            subTreeCosts.add(edgeCosts.get(i));
        }
        ArrayList<Double> prizePi = new ArrayList<Double>();
        for (int i : resultNodes) {
            prizePi.add(pi.get(i));
        }
        double totalPrizes = 0.0D;
        for (int i = 0; i < pi.size(); i++) {
            totalPrizes += pi.get(i);
        }
        HeadApprox pcsfHead = new HeadApprox();
        f = pcsfHead.new F(resultNodes, subTreeEdges, subTreeCosts, prizePi, totalPrizes);
    }


}
