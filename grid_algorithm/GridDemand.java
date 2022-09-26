package grid_algorithm;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import project_utils.Tuple;

/**
 * This class represents demands/flow divergence for grid graphs.
 * @author Jonas Schulz
 */
public class GridDemand implements Iterable<Tuple<Integer, Double>> {

    /**
     * Indices of relevant entries. 
     * "Irrelevant" currently means that the entry was never changed from the initial zero.
     * Relevant means not irrelevant.
     * Might be replaced or moved in future versions.
     */
    HashSet<Integer> entryIndices;
    /**
     * Entries of the demand/flow divergence b.
     * This stores the whole information of b.
     */
    double[] entries;
    /**
     * The graph that is related to this demand/flow divergence.
     */
    GridGraph g;
    
    /**
     * Standard Constructor.
     * Use this to construct the demand/flow divergence object, 
     * then set its entries via <code>set(Integer, Double)</code> or <code>set(int[], Double)</code>.
     * <i>Example call:</i> <code>new GridDemand(new GridGraph(4,4));</code>
     * @param g The graph related to the constructed demand/flow divergence.
     */
    public GridDemand(GridGraph g){
        this.g = g;
        entries = new double[g.getN()];
        entryIndices = new HashSet<>();
    }
    /**
     * Alternative Constructor.
     * Same as the standard constructor, but uses the graph provided by f.
     * Equivalent to <code>GridDemand(f.g)</code>.
     * @param f 
     */
    public GridDemand(GridFlow f){
        this.g = f.g;
        entries = new double[g.getN()];
        entryIndices = new HashSet<>();
    }
    /**
     * Deprecated Constructor.
     * @param grid_dimensions Dimensions of the grid graph this demand/flow divergence is related to.
     * @param n Number of nodes. Should equal &prod;<sub>i</sub>(<code>grid_dimensions[i]</code>).
     * @deprecated Functionalities moved to GridGraph class.
     */
    @Deprecated
    public GridDemand(int[] grid_dimensions, int n){
        this.g = new GridGraph(grid_dimensions);
        entries = new double[n];
        entryIndices = new HashSet<>();
    }
    
    /**
     * Sets the entry of b at index <code>at</code> to <code>value</code>.
     * @param at Index of the entry in vector notation (i.e. according to the enumeration scheme).
     * @param value Value that b<sub><code>at</code></sub> should be set to.
     */
    // currently not checking if contained to spare runtime
    public void set(Integer at, Double value){
        //entries.set(((int)at), value);
        entries[(int)at] = value;
        entryIndices.add(at);
    }
    /**
     * Sets the entry of b at index <code>at_index</code> to <code>value</code>.
     * Equivalent to <code>b.set(b.g.toIndex(at_index), value);
     * @param at_index Index in vector coordinates.
     * @param value Value that b<sub><code>at</code><sub> should be set to.
     */
    public void set(int[] at_index, Double value){
        set((int)GridGraph.toIndex(at_index, g.nodesPerDim), value);
    }
    /**
     * Adds <code>value</code> to b<sub><code>at</code></sub>.
     * @param at Index according to the enumeration scheme.
     * @param value Value to be added.
     */
    public void add(Integer at, Double value){
        if(entryIndices.contains(at)){
            //entries.set((int)at, entries.get((int)at) + value);
            entries[(int)at] += value;
        } else {
            entries[(int)at] += value;
            entryIndices.add(at);
        }
    }
    /**
     * Trivial Getter for b<sub><code>at</code></sub>.
     * @param at Index according to the enumeration scheme.
     * @return b<sub><code>at</code></sub>.
     */
    public double get(Integer at){
        return entries[(int)at];
    }
    
    /**
     * Trivial Getter for b<sub><code>at</code></sub>.
     * Equivalent to <code>b.get(b.g.toIndex(at));</code>.
     * @param at_index Index in vector coordinates.
     * @return b<sub><code>at</code></sub>.
     */
    public double get(int[] at_index){
        return entries[(int)GridGraph.toIndex(at_index, g.nodesPerDim)];
    }
    
    /**
     * Implements the <code>Iterator</code> interface.
     * @return The iterator that is used with <code>for (bi : b)</code>.
     */
    @Override
    public Iterator<Tuple<Integer, Double>> iterator() {
        return new Iterator<Tuple<Integer, Double>>() {
            Iterator<Integer> it = entryIndices.iterator();
            @Override
            public boolean hasNext() {
                return it.hasNext();
            }

            @Override
            public Tuple<Integer, Double> next() {
                Integer i = it.next();
                return new Tuple<>(i, entries[(int)i]);
            }
        };
    }
    
    /**
     * Calculates e<sup><i>lmax</i>(b)</sup>=&sum;<sub>i</sub> e<sup>b<sub>i</sub></sup>+e<sup>-b<sub>i</sub></sup>.
     * @return &sum;<sub>i</sub> e<sup>b<sub>i</sub></sup>+e<sup>-b<sub>i</sub></sup>.
     */
    public double lmax_exp(){
        double sum = 0;
        int non_zeros = 0;
        for(Integer e : entryIndices){
            double val = entries[(int)e];
            double exp = Math.exp(val);
            // choose following line iff Math.exp(x) is cheaper then 1/x
            //    sum += exp + Math.exp(-e.c);
            // else choose this line
            sum += exp + (1/exp);
            non_zeros++;
        }
        // as this is a sparse representation, but e^0 = 1, we have to add 
        // 2 times the amount of zero-entries
        sum += 2*(g.getN()-non_zeros);
        return sum;
    }
    /**
     * Calculates <i>lmax</i>(b)=<i>ln</i>(&sum;<sub>i</sub> e<sup>b<sub>i</sub></sup>+e<sup>-b<sub>i</sub></sup>).
     * @return <i>ln</i>(&sum;<sub>i</sub> e<sup>b<sub>i</sub></sup>+e<sup>-b<sub>i</sub></sup>).
     */
    public double lmax(){
        return Math.log(lmax_exp());
    }
    /**
     * Calculates the <i>l</i><sub>1</sub>-norm.
     * @return &Vert; b &Vert;<sub>1</sub> = &sum;<sub>i</sub> &vert; b<sub>i</sub> &vert;.
     */
    public double l1(){
        double sum = 0;
        for(Integer e : entryIndices){
            double val = Math.abs(entries[(int)e]);
            sum += val;
        }
        return sum;
    }
    /**
     * Calculates the <i>l</i><sub>&infin;</sub>-norm.
     * @return &Vert b &Vert;<sub>&infin;</sub> = <i>max</i><sub>i </sub>{&vert; b<sub>i</sub> &vert;}.
     */
    public double linf(){
        double max = 0;
        for(Integer e : entryIndices){
            double val = Math.abs(entries[(int)e]);
            if(val > max) max = val;
        }
        return max;
    }
    
    /**
     * Calculates the difference of two demands/flow divergences.
     * @param a The minuend.
     * @param b The subtrahend.
     * @return a - b.
     */
    public static GridDemand subtract(GridDemand a, GridDemand b){
        GridDemand res = new GridDemand(a.g);
        for(Tuple<Integer, Double> e : a){
            res.set((int)e.a, e.b);
        }
        for(Tuple<Integer, Double> e : b){
            res.add((int)e.a, -1*e.b);
        }
        return res;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected GridDemand clone(){
        GridDemand d = new GridDemand(g);
        for(Integer i : entryIndices){
            d.entryIndices.add(i);
            d.entries[i] = entries[i];
        }
        return d;
    }
    
    /**
     * Returns a copy of this demand/flow divergence, scaled by <code>s</code>.
     * @param s Scale factor.
     * @return <code>s</code> &#8729; b.
     */
    public GridDemand scale(double s) {
        GridDemand d = clone();
        for(Integer i : d.entryIndices){
            d.entries[i] *= s;
        }
        return d;
    }
    /**
     * Scales this demand/flow divergence by <code>s</code> (in-place).
     * @param s Scale factor.
     * @return <code>s</code> &#8729; b.
     */
    public GridDemand scale_inplace(double s) {
        for(Integer i : this.entryIndices){
            this.entries[i] *= s;
        }
        return this;
    }
    
    /**
     * Calculates the standard scalar product for two demands/flow divergences.
     * @param a First demand/flow divergence.
     * @param b Second demand/flow divergence.
     * @return &#x3008; <code>a</code> , <code>b</code> &#x3009; = <code>a</code><sup>T</sup> &#8729; <code>b</code>.
     */
    public static double scalar_prod(GridDemand a, GridDemand b){
        double sum = 0.;
        for(Integer i : a.entryIndices){
            sum += a.entries[i] * b.entries[i];
        }
        return sum;
    }
    
    /**
     * Calculates the flow given by the difference of the excess flows.
     * An edge (b<sub>i</sub>, b<sub>j</sub>) gets the value b<sub>j</sub>-b<sub>i</sub>.
     * @return B<sup>T</sup> &#8729; b.
     */
    public GridFlow toPotentialDiffEdgesFlow(){
        GridFlow f = new GridFlow(g);
        for(int i = 0; i < g.getN(); i++){
            int[][] neighbours = g.getNeighbours(g.toPosition(i));
            for(int k = 0; k < neighbours.length; k++){
                int other = g.toIndex(neighbours[k]);
                if(i < other){
                    double dif = this.get(i) - this.get(other);
                    f.set(i, other, -dif);
                }
            }
        }
        return f;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return Arrays.toString(this.entries);
    }
    
    /**
     * Returns TikZ-input for a visual representation of this demand/flow divergence.
     * Note that it currently <i>only</i> supports 2-dimensional grid graphs.
     * 1-dimensional graphs will result in a {@link java.lang.ArrayIndexOutOfBoundsException}.
     * Higher-dimensional graphs will throw all vertices with same coordinates in the first and second dimension
     * to one position, which results in TikZ only drawing the last vertex among those.
     * @return TikZ input for a visual representation of this demand/flow divergence.
     */
    public String tikz2D(){
        DecimalFormat df = new DecimalFormat("#.000");
        String s = "\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]\n";
        for(int i = 0; i < this.g.getN(); i++){
            double scale = 3;
            s += "\\node[roundnode, minimum size = 2cm] at ("+scale*g.toPosition(i)[0]+", "+scale*g.toPosition(i)[1]+") {\\textcolor{blue!28}{"+i+":} "+(get(i) == 0.0 ? "0" : "\\textcolor{red}{"+df.format(get(i))+"}")+"};\n";
        }
        s += "\\end{tikzpicture}";
        return s;
    }
    
}
