package grid_algorithm;

import java.util.ArrayList;
import jdk.jshell.spi.ExecutionControl;
import project_utils.*;

/**
 * This class represents a multidimensional grid graph with unit-capacity undirected edges.
 * @author Jonas Schulz
 */
public class GridGraph {
    /**
     * The vertices per dimension. For example, a 4-dimensional 2x4x5x3 grid would
     * have <code>nodesPerDim</code> = [2,4,5,3].
     */
    int[] nodesPerDim;
    
    /**
     * Standard constructor.
     * @param nodesPerDim Nodes per dimension.
     */
    public GridGraph(int... nodesPerDim){
        this.nodesPerDim = nodesPerDim;
    }
    
    /**
     * Returns the dimensionality of this grid graph.
     * Example: <code>(new GridGraph(2,3,2,3)).getDim()</code> returns 4.
     * @return Dimensionality of this grid graph.
     */
    public int getDim(){
        return nodesPerDim.length;
    }
    
    /**
     * Returns the number of vertices/nodes ("n")  of this grid graph.
     * @return n.
     */
    public int getN() {
        int prod = 1;
        for(int i : nodesPerDim){
            prod *= i;
        }
        return prod;
    }
    
    /**
     * Returns the number of edges ("m") of this grid graph.
     * The formula used to calculate m is: m = &sum;<sub>i</sub>(n/n<sub>i</sub>) &#8729; (n<sub>i</sub>-1).
     * @return m.
     */
    public int getM(){
        int sum = 0;
        int n = getN();
        for(int i : nodesPerDim){
            sum += (n/i)*(i-1);
        }
        return sum;
    }
    
    /**
     * Calculates the neighbour positions of a vertex.
     * @param at A vertex in vector grid coordinates.
     * @return The neighbours of <code>at</code>, as an array of vector grid coordinates.
     */
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
    
    /**
     * Calculates the cut capacity of a hypercube cut.
     * If <code>lowerBounds</code> is [l<sub>1</sub>,...,l<sub>d</sub>] and <code>higherBounds</code> is [h<sub>1</sub>,...,h<sub>d</sub>],
     * the represented hypercube contains all vertices [v<sub>1</sub>,...,v<sub>d</sub>] with l<sub>i</sub> &#8804; v<sub>i</sub> &#8804; h<sub>i</sub>
     * for all i in [d].
     * @param lowerBounds Lower Bounds of the hypercube.
     * @param higherBounds Higher Bounds of the hypercube.
     * @return Capacity of the hypercube cut.
     */
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
    
    /**
     * @deprecated Old implementation of an approximator functionality.
     */
    @Deprecated
    private double edgeScalar(int[] lowerBounds, int[] higherBounds, Vector v){
        int[] indices = indicesOfBoxNodes(lowerBounds, higherBounds);
        double s = 0;
        for(int i : indices){
            s += v.get(i);
        }
        s /= capHyperBox(lowerBounds, higherBounds);
        return s;
    }
    
    /**
     * Splits a hypercube into smaller hypercubes.<br>
     * The hypercube given in the input through its corners, <code>lowerBounds</code> and <code>higherBounds</code>,
     * is split at its center along each dimension that has more than one hyperplane.<br>
     * The result is an array of hypercube representations, where each hypercube representation
     * is an array with the two corners, each one in vector grid coordinates.<br>
     * For example, <code>split([2,0,2],[3,3,2])</code> will return<br>
     * <code>[<br>
     * &emsp;[&emsp;[2,0,2], [2,1,2]&emsp;],<br>
     * &emsp;[&emsp;[3,0,2], [3,1,2]&emsp;],<br>
     * &emsp;[&emsp;[2,2,2], [2,3,2]&emsp;],<br>
     * &emsp;[&emsp;[3,2,2], [3,3,2]&emsp;]<br>
     * ]</code>.
     * @param lowerBounds Lower bounds of the hypercube to be split.
     * @param higherBounds Higher bounds of the hypercube to be split.
     * @return An array with the hypercube splits of the input hypercube.
     */
    int[][][] split(int[] lowerBounds, int[] higherBounds){
        int dif = different(lowerBounds, higherBounds);
        int sc = countSplits(lowerBounds, higherBounds);
        int[][][] splits = new int[sc][2][lowerBounds.length];
        /**
         * The latter k iterates the actual dimensions, the latter i will denote the index of the currently calculated split.
         * To determine the split, we use the binary code of i.
         * As an example, let d = 5, and be the given hypercube 
         *     [[0,0,0,0,0], 
         *      [9,0,0,9,9]].
         *       ^ = = ^ ^
         *       0     3 4
         * This hypercube needs to be split along dimensions 0, 3 and 4.
         * Let us denote a split in the lower half of dimension 3 with 3l, the higher half with 3h (similar for other dimensions).
         * Example: [4l,3l,0h] for [[5,0,0,0,0],[9,0,0,4,4]].
         * i will iterate from 0 to 7 (countSplits([0,0,0,0,0],[9,0,0,9,9])==8).
         * We now try to do the following mapping:
         *   [4l,3l,0l]   [4l,3l,0h]   [4l,3h,0l]   [4l,3h,0h]   [4h,3l,0l]   [4h,3l,0h]   [4h,3h,0l]   [4h,3h,0h]   (hypercubes)
         *     0  0  0      0  0  1      0  1  0      0  1  1      1  0  0      1  0  1      1  1  0      1  1  1    (i - binary)
         *        0            1            2            3            4            5            6            7       (i - decimal)
         * We introduce the map dim_indices to store the dimensions with a split along them, and its inverse dim_ind_inverse.
         * For example, dim_indices == [0, 3, 4] and dim_ind_inverse == [0, -1, -1, 1, 2].
         */
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
            for(int k = 0; k < lowerBounds.length; k++){
                if(lowerBounds[k] < higherBounds[k]){
                    /** 
                     * We now do the look-up with the binary code of i.
                     * With ib as binary code of i, we should currently handle the lower split iff 
                     * the dimension k is denoted with "l" (or "h" for upper split) in the description above.
                     * We get this information via the dim_ind_inverse[k]-th bit of ib.
                     * This gives us the following look-up.
                     */
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
    
    /**
     * Counts the splits for <code>split(lowerBounds, higherBounds)</code>.
     * Example: <code>countSplits([0,0,0,0,0],[9,0,0,9,9])</code> will return 2<sup>3</sup> = 8.
     * @param lowerBounds Lower bounds of the represented hypercube.
     * @param higherBounds Higher Bounds of the represented hypercube.
     * @return Number of splits that <code>split(lowerBounds, higherBounds)</code> will use.
     */
    int countSplits(int[] lowerBounds, int[] higherBounds){
        int different = different(lowerBounds, higherBounds);
        return 1 << different;
    }
    
    /**
     * Counts the number of different indices of <code>a</code> and <code>b</code>.
     * Only use this method if both arrays have same length.
     * Example: <code>different([0,0,0,0,0],[9,0,0,9,9])</code> will return 3.
     * @param a First array for comparison.
     * @param b Second array for comparison.
     * @return Number of different indices.
     */
    int different(int[] a, int[] b){
        int dif = 0;
        for(int i = 0; i < a.length; i++){
            if(a[i] != b[i]) dif++;
        }
        return dif;
    }
    
    /**
     * Recursively calculates the indices according to the enumeration scheme for all vertices inside the given hypercube.
     * It is coherent with the enumeration scheme from <code>toIndex</code>.
     * @param lowerBounds Lower bounds of the represented hypercube.
     * @param higherBounds Higher bounds of the represented hypercube.
     * @return Array with the indices of all vertices inside the hypercube, according to the enumeration scheme from <code>toIndex</code>.
     */
    public int[] indicesOfBoxNodes(int[] lowerBounds, int[] higherBounds){
        int iterationLevel = 0;
        while(iterationLevel < nodesPerDim.length && higherBounds[iterationLevel] == lowerBounds[iterationLevel]){
            iterationLevel++;
        }
        /**
         * Determines the current recursive iteration level.
         * The first index where lower and higher bounds do not differ is the dimension along which we
         * recursively call this method. For example:
         *      indicesOfBoxNodes([2,0,2,2,4,1], [2,0,2,4,7,1])
         *                               ^              ^
         *   will make recursive calls to 
         *      indicesOfBoxNodes([2,0,2,2,4,1], [2,0,2,2,7,1])
         *      indicesOfBoxNodes([2,0,2,3,4,1], [2,0,2,3,7,1])
         *      indicesOfBoxNodes([2,0,2,4,4,1], [2,0,2,4,7,1])
         * iterationLevel == nodesPerDim.length  IFF  lowerBounds == higherBounds
         */
        if(iterationLevel >= nodesPerDim.length) return new int[]{toIndex(lowerBounds)};
        /**
         * iterate recursively the current level and add the results to "indices"
         */
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
    
    /**
     * Calculates the volume of a hypercube.
     * The volume of a hypercube is the number of vertices contained in it.
     * @param lowerIndices Lower bounds of the hypercube.
     * @param higherIndices Higher bounds of the hypercube.
     * @return The volume of the given hypercube.
     */
    public int volume(int[] lowerIndices, int[] higherIndices){
        int prod = 1;
        for(int i = 0; i < lowerIndices.length; i++){
            prod *= higherIndices[i] - lowerIndices[i] + 1;
        }
        return prod;
    }
    
    /**
     * Defines the enumeration of coordinate vectors inside this grid graph.<br>
     * This function implements a bijection from [n<sub>1</sub>]&times;...&times;[n<sub>d</sub>]
     * (grid graph coordinate vector space) to [n] = [n<sub>1</sub>&#8729;...&#8729;n<sub>d</sub>] (natural numbers).<br>
     * Note that [x] is adjusted in the code to [x-1]<sub>0</sub>, as usual for implementations of mathematical objects like vectors.<br>
     * <br>
     * Example: Let G be a grid graph with dimensions <code>[2,3,3]</code>.<br>
     * This function will yield the following bijection for G:<br>
     * <code>
     * &emsp; toIndex[0,0,0] ==  0<br>
     * &emsp; toIndex[0,0,1] ==  1<br>
     * &emsp; toIndex[0,0,2] ==  2<br>
     * &emsp; toIndex[0,1,0] ==  3<br>
     * &emsp; toIndex[0,1,1] ==  4<br>
     * &emsp; toIndex[0,1,2] ==  5<br>
     * &emsp; toIndex[0,2,0] ==  6<br>
     * &emsp; toIndex[0,2,1] ==  7<br>
     * &emsp; toIndex[0,2,2] ==  8<br>
     * &emsp; toIndex[1,0,0] ==  9<br>
     * &emsp; toIndex[1,0,1] == 10<br>
     * &emsp; toIndex[1,0,2] == 11<br>
     * &emsp; toIndex[1,1,0] == 12<br>
     * &emsp; toIndex[1,1,1] == 13<br>
     * &emsp; toIndex[1,1,2] == 14<br>
     * &emsp; toIndex[1,2,0] == 15<br>
     * &emsp; toIndex[1,2,1] == 16<br>
     * &emsp; toIndex[1,2,2] == 17<br>
     * </code>
     * @param position Grid graph coordinate vector.
     * @return Enumerated index of <code>position</code>.
     */
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
    
    /**
     * Same as <code>toIndex(int[])</code>, but enumerates a general finite set over &#8469;<sup>d</sup>.
     * @param position Vector in [s<sub>1</sub>]&times;...&times;[s<sub>d</sub>], where s<sub>i</sub> is <code>sizePerDim[i]</code>.
     * @param sizePerDim [|S<sub>1</sub>|,...,|S<sub>d</sub>|], where the function domain is S<sub>1</sub>&times;...&times;S<sub>d</sub> with S<sub>i</sub>=[s<sub>i</sub>].
     * @return Enumerated index of <code>position</code>.
     */
    public static int toIndex(int[] position, int[] sizePerDim){
        assert position.length == sizePerDim.length;
        int index = 0;
        int fac = 1;
        for(int i = position.length - 1; i >= 0; i--){
            index += position[i]*fac;
            fac *= sizePerDim[i];
        }
        return index;
    }
    
    /**
     * Inverse function of <code>toIndex</code>.<br>
     * It holds that <br>
     * <code>&emsp;toIndex(toPosition(i)) == i</code> and<br>
     * <code>&emsp;toPosition(toIndex(v)) == v</code>.
     * @param index Index in [n] (or code notation: [n-1]<sub>0</sub>).
     * @return Grid graph coordinate vector v that satisfies <code>toIndex(v) == index</code>.
     */
    // Tested
    public int[] toPosition(int index){
        int[] pos = new int[nodesPerDim.length];
        for(int i = 0; i < nodesPerDim.length; i++){
            int j = nodesPerDim.length - i - 1;
            int current_div = nodesPerDim[j];
            pos[j] = index - (index / current_div)*current_div;
            index = (index - pos[j]) / current_div;
        }
        return pos;
    }
    
    /**
     * Same as <code>toPosition</code>, but uses Java's built-in modulo symbol.
     */
    public int[] toPositionV2(int index){
        int[] pos = new int[nodesPerDim.length];
        for(int j = nodesPerDim.length - 1; j >= 0; j--){
            int current_div = nodesPerDim[j];
            pos[j] = index % current_div;
            index = (index - pos[j]) / current_div;
        }
        return pos;
    }
    
    /**
     * Checks if all given coordinate vectors are contained in the domain of this grid graph.
     * @param bounds Array containing the coordinate vectors to be checked.
     * @return <code>true</code> iff <b>all</b> given coordinate vectors are contained in the domain of this grid graph.
     */
    public boolean inbound(int[][] bounds){
        for(int i = 0; i < nodesPerDim.length; i++){
            for(int j = 0; j < bounds.length; j++){
                if(bounds[j][i] >= nodesPerDim[i]) return false;
            }
        }
        return true;
    }
    
    /**
     * Returns TikZ-input for a visual representation of this grid graph.
     * Note that it currently <i>only</i> supports 2-dimensional grid graphs.
     * 1-dimensional graphs will result in a {@link java.lang.ArrayIndexOutOfBoundsException}.
     * Higher-dimensional graphs will throw all vertices with same coordinates in the first and second dimension
     * to one position, which results in TikZ only drawing the last vertex among those.
     * @return TikZ input for a visual representation of this grid graph.
     */
    public String tikz2D(){
        String s = "";
        //String s = "\\begin{tikzpicture}\n";
        //String s = "\\node (anchor) {};\n";
        for(int i = 0; i < getN(); i++){
            double scale = 2;
            s += "\\node[roundnode, minimum size = 1cm] at ("+scale*toPosition(i)[0]+", "+scale*toPosition(i)[1]+") ("+i+") {\\textcolor{black}{"+i+"}};\n";
            int[][] neighbours = getNeighbours(toPosition(i));
            for(int[] n : neighbours){
                int indexTo = toIndex(n);
                if(indexTo < i){
                    s += "\\draw ("+i+") -- ("+indexTo+");\n";
                }
            }
        }
        //s += "\\end{tikzpicture}";
        return s;
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
