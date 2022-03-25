/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grid_algorithm;

import java.util.ArrayList;
import jdk.jshell.spi.ExecutionControl;
import project_utils.*;

/**
 *
 * @author kroka
 */
public class GridGraph /*implements RoutingGraph*/ {
    int[] nodesPerDim;
    public GridGraph(int... nodesPerDim){
        this.nodesPerDim = nodesPerDim;
    }
    public int getDim(){
        return nodesPerDim.length;
    }
    public int getN() {
        int prod = 1;
        for(int i : nodesPerDim){
            prod *= i;
        }
        return prod;
    }
    // todo: make more efficient 1d-index variant, which might also reduce
    // loop complexities to theta(d) instead theta(d^2) through index calculations
    // Tested
    public int[][] getNeighbours(int[] at) {
        int size = 2*at.length;
        for(int i = 0; i < at.length; i++){
            if(at[i] == 0)
                size--;
            if(at[i] == nodesPerDim[i] - 1)
                size--;
        }
        int[][] neighbours = new int[size][];
        for(int i = 0, k = 0; i < at.length; i++){
            if(at[i] != 0){
                neighbours[k] = new int[at.length];
                for(int j = 0; j < at.length; j++){
                    if(i == j){
                        neighbours[k][j] = at[j] - 1;
                    } else {
                        neighbours[k][j] = at[j];
                    }
                }
                k++;
            }
            if(at[i] != nodesPerDim[i] - 1){
                neighbours[k] = new int[at.length];
                for(int j = 0; j < at.length; j++){
                    if(i == j){
                        neighbours[k][j] = at[j] + 1;
                    } else {
                        neighbours[k][j] = at[j];
                    }
                }
                k++;
            }
        }
        return neighbours;
    }

//    @Deprecated
//    @Override
//    public Vector approximateCongestion(Vector v) {
//        assert v.getN() == this.getN();
//        ArrayList<Double> edgeScalars = new ArrayList<>();
//        int[] lower = new int[nodesPerDim.length];
//        int[] higher = new int[nodesPerDim.length];
//        for(int i = 0; i < nodesPerDim.length; i++){
//            lower[i] = 0;
//            higher[i] = nodesPerDim[i] - 1;
//        }
//        addEdgeScalarToListRecursive(new int[][]{lower, higher}, edgeScalars, v, false);
//        double[] ds = new double[edgeScalars.size()];
//        for(int i = 0; i < edgeScalars.size(); i++){
//            ds[i] = edgeScalars.get(i);
//        }
//        return new ArrayVector(ds);
//    }
//    @Deprecated
//    private void addEdgeScalarToListRecursive(int[][] bounds, ArrayList<Double> list, Vector v, boolean add){
//        if(add)
//            list.add(edgeScalar(bounds[0], bounds[1], v));
//        if(different(bounds[0], bounds[1]) != 0){
//            int[][][] splits = split(bounds[0], bounds[1]);
//            for(int[][] ia : splits){
//                addEdgeScalarToListRecursive(new int[][]{ia[0], ia[1]}, list, v, true);
//            }
//        }
//    }
    
    public int capHyperBox(int[] lowerBounds, int[] higherBounds){
        assert lowerBounds.length == higherBounds.length && lowerBounds.length == nodesPerDim.length;
        int total = 0;
        int prodAll = 1;
        for(int i = 0; i < lowerBounds.length; i++){
            prodAll *= higherBounds[i] - lowerBounds[i] + 1;
        }
        for (int i = 0; i < nodesPerDim.length; i++) {
            int fac = 0;
            if (lowerBounds[i] != 0) fac++;
            if (higherBounds[i] != nodesPerDim[i]-1) fac++;
            if (fac != 0){
                total += fac*prodAll/(higherBounds[i]-lowerBounds[i]+1);
            }
        }
        return total;
    }
    private double edgeScalar(int[] lowerBounds, int[] higherBounds, Vector v){
        int[] indices = indicesOfBoxNodes(lowerBounds, higherBounds);
        double s = 0;
        for(int i : indices){
            s += v.get(i);
        }
        s /= capHyperBox(lowerBounds, higherBounds);
        return s;
    }
    int[][][] split(int[] lowerBounds, int[] higherBounds){
        int dif = different(lowerBounds, higherBounds);
        int sc = countSplits(lowerBounds, higherBounds);
        int[][][] splits = new int[sc][2][lowerBounds.length];
        int[] dim_indices = new int[dif];
        int[] dim_ind_inverse = new int[lowerBounds.length];
        for(int i = 0, k = 0; i < lowerBounds.length; i++){
            if(lowerBounds[i] < higherBounds[i]){
                dim_ind_inverse[i] = k;
                dim_indices[k++] = i;
            } else {
                dim_ind_inverse[i] = -1;
            }
        }
        for(int i = 0; i < sc; i++){
            // each bit is representing one dimension (that is different in low/high index)
            // 0 can be assigned the lower half, 1 the bigger half
            for(int k = 0; k < lowerBounds.length; k++){
                if(lowerBounds[k] < higherBounds[k]){
                    assert dim_ind_inverse[k] != -1;
                    boolean low = (i & (1 << dim_ind_inverse[k])) == 0;
                    int med = (lowerBounds[k] + higherBounds[k])/2;
                    splits[i][0][k] = low ? lowerBounds[k] : (med + 1);
                    splits[i][1][k] = low ? med : higherBounds[k];
                } else {
                    // lower/higher bounds are the same
                    splits[i][0][k] = lowerBounds[k];
                    splits[i][1][k] = lowerBounds[k];
                }
            }
        }
        return splits;
    }
    int countSplits(int[] lowerBounds, int[] higherBounds){
        int different = different(lowerBounds, higherBounds);
        return 1 << different;
    }
    int different(int[] a, int[] b){
        int dif = 0;
        for(int i = 0; i < a.length; i++){
            if(a[i] != b[i]) dif++;
        }
        return dif;
    }
    public int[] indicesOfBoxNodes(int[] lowerBounds, int[] higherBounds){
        int iterationLevel = 0;
        while(iterationLevel < nodesPerDim.length && higherBounds[iterationLevel] == lowerBounds[iterationLevel]){
            iterationLevel++;
        }
        // iterationLevel = nodesPerDim.length <-> lowerBounds == higherBounds
        if(iterationLevel >= nodesPerDim.length) return new int[]{toIndex(lowerBounds)};
        // iterate recursively the current level and add the results to indices
        ArrayList<Integer> indices = new ArrayList<>();
        int[] currentRecursionLower = new int[lowerBounds.length];
        int[] currentRecursionHigher = new int[higherBounds.length];
        for(int i = 0; i < lowerBounds.length; i++){
            currentRecursionLower[i] = lowerBounds[i];
            currentRecursionHigher[i] = higherBounds[i];
        }
        for(int i = lowerBounds[iterationLevel]; i <= higherBounds[iterationLevel]; i++){
            currentRecursionLower[iterationLevel] = i;
            currentRecursionHigher[iterationLevel] = i;
            int[] iterationLevelIndices = indicesOfBoxNodes(currentRecursionLower, currentRecursionHigher);
            for(int k : iterationLevelIndices) indices.add(k);
        }
        int[] result = new int[indices.size()];
        for(int k = 0; k < indices.size(); k++){
            result[k] = indices.get(k);
        }
        return result;
    }
    public int volume(int[] lowerIndices, int[] higherIndices){
        int prod = 1;
        for(int i = 0; i < lowerIndices.length; i++){
            prod *= higherIndices[i] - lowerIndices[i] + 1;
        }
        return prod;
    }
    // Tested
    public int toIndex(int[] position){
        assert position.length == nodesPerDim.length;
        int index = 0;
        int fac = 1;
        for(int i = position.length - 1; i >= 0; i--){
            index += position[i]*fac;
            fac *= nodesPerDim[i];
        }
        return index;
    }
    public static int toIndex(int[] position, int[] nodesPerDim){
        assert position.length == nodesPerDim.length;
        int index = 0;
        int fac = 1;
        for(int i = position.length - 1; i >= 0; i--){
            index += position[i]*fac;
            fac *= nodesPerDim[i];
        }
        return index;
    }
    // Tested
    public int[] toPosition(int index){
        int[] pos = new int[nodesPerDim.length];
        // int prev_index = index;
        // int current_div = 1;
        for(int i = 0; i < nodesPerDim.length; i++){
            int j = nodesPerDim.length - i - 1;
            // current_div *= nodesPerDim[j];
            int current_div = nodesPerDim[j];
            pos[j] = index - (index / current_div)*current_div;
            index = (index - pos[j]) / current_div;
            
            // prev_index = index;
            // index = index / current_div;
            // pos[j] = prev_index - index*current_div;
            // index -= pos[j];
        }
        return pos;
    }
    
    public boolean inbound(int[][] bounds){
        for(int i = 0; i < nodesPerDim.length; i++){
            for(int j = 0; j < bounds.length; j++){
                if(bounds[j][i] >= nodesPerDim[i]) return false;
            }
        }
        return true;
    }
//    public double getAlpha(){
//        return getN();
//    }
    
//    public HashTree getApproximatorAsTree(){
//        
//    }
    
//    @Deprecated
//    public Vector demand(Vector flow){
//        return flow;
//    }
//    @Deprecated
//    public double potential(Vector flow, Vector demands){
//        double pot_graph = flow.symm_softmax();
//        Vector dem = demand(flow);
//        Vector residual = demands.add(dem.scalarMultiply(-1));
//        Vector vec_tree = approximateCongestion(residual);
//        return pot_graph + (vec_tree.scalarMultiply(2/* /*getAlpha()*/)).symm_softmax();
//    }
}
