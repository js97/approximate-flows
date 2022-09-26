package grid_algorithm;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map.Entry;
import project_utils.Triple;
import project_utils.Tuple;

/**
 * This class represents flows for grid graphs.
 * @author Jonas Schulz
 */
public class GridFlow {
    /**
     * The grid graph related to this flow.
     */
    GridGraph g;
    /**
     * A {@link java.util.HashMap} storing the information of this flow.
     * Might be changed to a {@link java.lang.Long} &rarr; {@link java.lang.Double} map for bigger graph support.
     * Might also be replaced by an array-based implementation, once an enumeration scheme (and its inverse) is established.
     */
    HashMap<Integer, Double> entries;
    /**
     * Stores the number of edges.
     * With the current implementation, the edge count has to be calculated. 
     * To improve efficiency, it is calculated once at construction (or assigned) and can then be used.
     * Will be replaced as soon as the entries are stored in an array with length m.
     */
    int m;
    
    /**
     * Standard constructor.
     * Calculates m and initializes the fields.
     * @param g The graph related to this flow.
     */
    public GridFlow(GridGraph g){
        this.g = g;
        calculateM();
        entries = new HashMap<>(m);
    }
    
    /**
     * Alternative "copy"/"mirror" constructor.
     * The constructed flow object will be the same as the input, but only with 0 at each edge iff clone is false.
     * Iff clone is true, the constructed flow object is a copy of f.
     * @param f Flow object to clone/mirror.
     * @param clone <ul><li><code>true</code> - construct a clone of <code>f</code>.</li><li><code>false</code> - construct a flow with same structure as <code>f</code>.</li></ul>
     */
    GridFlow(GridFlow f, boolean clone){
        this.entries = new HashMap<>();
        this.g = f.g;
        this.m = f.m;
        if(clone){
            for(Entry<Integer, Double> t : f.entries.entrySet()){
                this.entries.put(t.getKey(), t.getValue());
            }
        }
    }
    
    /**
     * Calculates m.
     * The formula used to calculate m is: m = &sum;<sub>i</sub>(n/n<sub>i</sub>) &#8729; (n<sub>i</sub>-1).
     * @return m.
     */
    public int calculateM(){
        int sum = 0;
        for(int i = 0; i < g.nodesPerDim.length; i++){
            sum += (g.getN()/g.nodesPerDim[i])*(g.nodesPerDim[i]-1);
        }
        // 2*sum or sum ? TODO
        this.m = sum;
        return this.m;
    }

    /**
     * Standard Getter for m (edge count).
     * @return m.
     */
    public int getM() {
        return m;
    }

    /**
     * Getter for n (vertex count).
     * @return n.
     */
    public int getN() {
        return g.getN();
    }
    
    /**
     * Set f<sub>e</sub>= <code>value</code>, where e = (<code>from,to</code>).
     * @param from First vertex of the edge.
     * @param to Second vertex of the edge.
     * @param value Value to be assigned to the edge.
     */
    public void set(Integer from, Integer to, Double value){
        Triple<Integer, Integer, Double> e = (from > to) ? new Triple<>(to, from, -value) : new Triple<>(from, to, value);
        entries.put(this.hashEntry(new Tuple<>(e.a, e.b)), e.c);
    }
    /**
     * Set f<sub>e</sub>= <code>value</code>, where e = (<code>from,to</code>).
     * Same as the other <code>set</code> method, but with grid coordinates for the vertices.
     * @param from First vertex of the edge.
     * @param to Second vertex of the edge.
     * @param value Value to be assigned to the edge.
     */
    public void set(int[] from, int[] to, Double value){
        int from_int = g.toIndex(from);
        int to_int = g.toIndex(to);
        this.set(from_int, to_int, value);
    }
    /**
     * Add <code>value</code> to f<sub>e</sub>, where e = (<code>from,to</code>).
     * @param from First vertex of the edge.
     * @param to Second vertex of the edge.
     * @param value Value to be added to the edge.
     */
    public void add(Integer from, Integer to, double value){
        if(from > to)
            this.add(this.hashEntry(new Tuple<>(to, from)), -value);
        else
            this.add(this.hashEntry(new Tuple<>(from, to)), value);
    }
    /**
     * Set f<sub>e</sub>= <code>value</code>, where <code>combinedIndex</code> is the identifier of e.
     * @param combinedIndex Identifier for the edge.
     * @param value Value to be assigned to the edge.
     */
    void set(int combinedIndex, double value){
        entries.put(combinedIndex, value);
    }
    /**
     * Add <code>value</code> to f<sub>e</sub>, where <code>combinedIndex</code> is the identifier of e.
     * @param combinedIndex Identifier for the edge.
     * @param value Value to be added to the edge.
     */
    void add(int combinedIndex, double value){
        entries.put(combinedIndex, (entries.containsKey(combinedIndex) ? entries.get(combinedIndex) : 0)+ value);
    }
    /**
     * Get f<sub>e</sub>, where e = (<code>from,to</code>).
     * @param from First vertex of the edge.
     * @param to Second vertex of the edge.
     * @return f<sub>e</sub>.
     */
    public double get(Integer from, Integer to){
        int combined = this.hashEntry((from > to) ? (new Tuple<>(to, from)) : (new Tuple<>(from, to)));
        return (entries.containsKey(combined)) ? ((from > to ? -1 : 1)*entries.get(combined)) : 0.;
    }
    /**
     * Get f<sub>e</sub>, where e = (<code>from,to</code>).
     * Same as the other <code>get</code> method, but with grid coordinates for the vertices.
     * @param from First vertex of the edge.
     * @param to Second vertex of the edge.
     * @return f<sub>e</sub>.
     */
    public double get(int[] from, int[] to){
        return get(this.g.toIndex(from), this.g.toIndex(to));
    }
    
    /**
     * Calculates <i>B &#8729; f</i> with:
     * <ul>
     * <li> <i>B</i> - the incidence matrix,</li>
     * <li> <i>f</i> - the flow.</li>
     * </ul>
     * @return <i>B &#8729; f</i>
     */
    public GridDemand calculateExcessFlows(){
        GridDemand b = new GridDemand(this);
        for(Entry<Integer, Double> e : entries.entrySet()){
            Tuple<Integer, Integer> edge = this.hashToTuple(e.getKey());
            b.add(edge.a, -e.getValue());
            b.add(edge.b, e.getValue());
        }
        return b;
    }
    
    /**
     * Calculates e<sup><i>lmax</i>(f)</sup>=&sum;<sub>i</sub> e<sup>f<sub>i</sub></sup>+e<sup>-f<sub>i</sub></sup>.
     * If some absolute value in <i>f</i> is greater than or equal to 709.7827128933841, this method is guaranteed to return Infinity.
     * To prevent this, use the shifted version instead.
     * Also, 709.782712893384 can be calculated in exp(x), but adding more values of this size also results in evaluating Infinity.
     * It should hold that m &#8729; exp(f<sub>i</sub>) won't exceed exp(709.782712893384).
     * Thus, if all entry absolutes in this flow are below (709.782712893384-ln(m)), this method should still be able to perform the calculation.
     * @return &sum;<sub>i</sub> e<sup>f<sub>i</sub></sup>+e<sup>-f<sub>i</sub></sup>.
     */
    public double lmax_exp(){
        double sum = 0;
        int non_zeros = 0;
        for(Entry<Integer, Double> e : entries.entrySet()){
            double exp = Math.exp(e.getValue());
            // choose following line iff Math.exp(-x) is cheaper then 1/x
            // sum += exp + Math.exp(-e.c);
            // else choose this line
            sum += exp + (1/exp);
            non_zeros++;
        }
        // as this is a sparse representation, but e^0 = 1, we have to add 
        // 2 times the amount of zero-entries
        sum += 2*(m-non_zeros);
        return sum;
    }
    /**
     * Calculates e<sup><i>lmax</i>(f)</sup>=&sum;<sub>i</sub> e<sup>f<sub>i</sub>+<code>shift</code></sup>+e<sup>-f<sub>i</sub>+<code>shift</code></sup>.
     * Use this method to avoid numeric errors due to high values in <i>f</i>. 
     * The shift is directly added to all exponents, so use a negative value for <code>shift</code> to decrease the exponents.
     * @return &sum;<sub>i</sub> e<sup>f<sub>i</sub>+<code>shift</code></sup>+e<sup>-f<sub>i</sub>+<code>shift</code></sup>.
     */
    public double lmax_exp_shifted(double shift){
        double sum_exps = 0;
        int nonzeros = 0;
        for(Entry<Integer, Double> e : entries.entrySet()){
            double exp = Math.exp(e.getValue() + shift) + Math.exp((-e.getValue()) + shift);
            sum_exps += exp;
            nonzeros++;
        }
        double exp_shift = Math.exp(shift);
        sum_exps += 2 * exp_shift * (m - nonzeros);
        return sum_exps;
    }
    /**
     * Calculates <i>lmax(f)</i>.
     * This version does not use shifts.
     * @return <i>lmax(f)</i>.
     */
    public double lmax(){
        // might use shifted version with shift = 0
        return Math.log(lmax_exp());
    }
    /**
     * Calculates <i>lmax(f)</i>.
     * This version shifts the exponents according to the parameter.
     * The mathematical calculation yields the same as the variant without shifts,
     * but this version prevents overflows.
     * @param shift The shift to be <b>added</b> to all exponents.
     * @return <i>lmax(f)</i>.
     */
    public double lmax_shifted(double shift){
        return Math.log(lmax_exp_shifted(shift)) - shift;
    }
    /**
     * Calculates <i>lmax(f)</i>.
     * Same as <code>lmax_shifted(shift)</code>, but uses the default shift.
     * @return <i>lmax(f)</i>.
     */
    public double lmax_shifted(){
        double shift = getDefaultShift();
        return lmax_shifted(shift);
    }
    /**
     * Calculates &nabla; <i>lmax(f)</i>.
     * Use the shifted variant if your exponents might cause numeric overflow.
     * @return &nabla; <i>lmax(f)</i>.
     */
    public GridFlow gradient_lmax(){
        double sum_exps = lmax_exp();
        GridFlow f = new GridFlow(this, false);
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            double exp = Math.exp(t.getValue());
            f.set(t.getKey(), (exp-(1./exp))/sum_exps);
        }
        // where there is i with no entry, xi = 0, hence e^xi - e^(-xi) = 0, thus 
        // also the gradient is 0 at i
        return f;
    }
    /**
     * Calculates &nabla; <i>lmax(f)</i>.
     * This version uses shifts, preventing numeric overflow to Infinity.
     * Use a negative value for <code>shift</code> to reduce exponents.
     * @return &nabla; <i>lmax(f)</i>.
     */
    public GridFlow gradient_lmax_shifted(double shift){
        double sum_exps = lmax_exp_shifted(shift);
        GridFlow f = new GridFlow(this, false);
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            double exp1 = Math.exp(t.getValue() + shift);
            double exp2 = Math.exp(-t.getValue() + shift);
            f.set(t.getKey(), (exp1 - exp2) / sum_exps);
        }
        return f;
    }
    /**
     * Calculates the default shift.
     * The default shift is the negative of the highest of all exponents, which is equivalent to
     * <i>-1 &#8729; &Vert; f &Vert;<sub>&infin;</sub></i>.
     * This shifts the exponents s.t. the maximum shifted exponent is 0.
     * @return the default shift of <i>-1 &#8729; &Vert; f &Vert;<sub>&infin;</sub></i>.
     */
    double getDefaultShift(){
        return -1*linf();
    }
    /**
     * Calculates &nabla; <i>lmax(f)</i>.
     * This version uses shifts, preventing numeric overflow to Infinity.
     * All shifts are set to the default shift.
     * @return &nabla; <i>lmax(f)</i>.
     */
    public GridFlow gradient_lmax_shifted(){
        return gradient_lmax_shifted(getDefaultShift());
    }
    /**
     * Calculates the <i>l</i><sub>1</sub>-norm.
     * @return &Vert; f &Vert;<sub>1</sub> = &sum;<sub>i</sub> &vert; f<sub>i</sub> &vert;.
     */
    public double l1() {
        double curr = 0;
        for(Entry<Integer, Double> e : entries.entrySet()){
            curr += Math.abs(e.getValue());
        }
        return curr;
    }
    /**
     * Calculates the <i>l</i><sub>&infin;</sub>-norm.
     * @return &Vert f &Vert;<sub>&infin;</sub> = <i>max</i><sub>i </sub>{&vert; f<sub>i</sub> &vert;}.
     */
    public double linf(){
        double max = 0;
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            double v = Math.abs(t.getValue());
            if(v > max) max = v;
        }
        return max;
    }

    /**
     * Returns a copy of this flow, scaled by <code>d</code>.
     * @param d Scale factor.
     * @return <code>d</code> &#8729; f.
     */
    public GridFlow scale(double d) {
        GridFlow cl = new GridFlow(this, true);
        for(Entry<Integer, Double> t : cl.entries.entrySet()){
            t.setValue(t.getValue() * d);
        }
        return cl;
    }
    /**
     * Scales this flow by <code>d</code> (in-place).
     * @param d Scale factor.
     * @return <code>d</code> &#8729; f.
     */
    public GridFlow scale_inplace(double d) {
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            t.setValue(t.getValue() * d);
        }
        return this;
    }
    /**
     * Adds two flows.
     * @param a First summand.
     * @param b Second summand.
     * @return <code>a</code> + <code>b</code>.
     */
    public static GridFlow add(GridFlow a, GridFlow b) {
        GridFlow cl = new GridFlow(a, true);
        for(Entry<Integer, Double> t : b.entries.entrySet()){
            cl.add(t.getKey(), t.getValue());
        }
        return cl;
    }
    
    /**
     * Returns the identifier for the edge <code>t</code>.
     * Might be extended to a <code>long</code> identifier in the future, or
     * even be replaced when an enumeration scheme for the edges has been established.
     * @param t The edge represented as tuple of vertex identifiers.
     * @return The identifier for this edge.
     */
    int hashEntry(Tuple<Integer, Integer> t){
        return (((int)t.a) << 16) + t.b;
    }
    /**
     * Inverse function to <code>hashEntry</code>.
     * @param hash Edge identifier.
     * @return Tuple of vertex identifiers.
     */
    Tuple<Integer, Integer> hashToTuple(int hash){
        int second = hash & 0x0000_ffff;
        int first = (hash >> 16) & 0x0000_ffff;
        return new Tuple<>(first, second);
    }

    /**
     * {@inheritDoc }
     */
    @Override
    public String toString() { 
        String s = "[\n";
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            Tuple<Integer, Integer> i = this.hashToTuple(t.getKey());
            s += "\t"+i.toString()+" :   "+t.getValue()+"\n";
        }
        return s + "]";
    }
    
    /**
     * Returns TikZ-input for a visual representation of this flow.
     * Note that it currently <i>only</i> supports 2-dimensional grid graphs.
     * 1-dimensional graphs will result in a {@link java.lang.ArrayIndexOutOfBoundsException}.
     * Higher-dimensional graphs will throw all vertices with same coordinates in the first and second dimension
     * to one position, which results in TikZ only drawing the last vertex among those.
     * @return TikZ input for a visual representation of this flow.
     */
    public String tikz2D(){
        String s = "\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]\n";
        DecimalFormat df = new DecimalFormat("#.000");
        for(int i = 0; i < getN(); i++){
            double scale = 2;
            s += "\\node[roundnode, minimum size = 1cm] at ("+scale*g.toPosition(i)[0]+", "+scale*g.toPosition(i)[1]+") ("+i+") {\\textcolor{black}{"+i+"}};\n";
        }
        for (int i = 0; i < getN(); i++) {
            int[][] neighbours = g.getNeighbours(g.toPosition(i));
            for(int[] n : neighbours){
                int indexTo = g.toIndex(n);
                if(indexTo > i){
                    String number = (this.get(i, indexTo) == 0.0) ? "\\textcolor{blue!28}{0}" : df.format(this.get(i, indexTo));
                    if(indexTo == i + 1)
                        s += "\\draw[thick, ->] ("+i+") -- ("+indexTo+") node[midway,right] {"+number+"};\n";
                    else
                        s += "\\draw[thick, ->] ("+i+") -- ("+indexTo+") node[midway,above] {"+number+"};\n";
                }
            }
            
        }
        s += "\\end{tikzpicture}";
        return s;
    }
}
