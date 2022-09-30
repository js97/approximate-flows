package grid_algorithm;

import java.awt.Color;
import java.text.DecimalFormat;
import java.util.Arrays;

import org.scilab.forge.jlatexmath.TeXConstants;
import org.scilab.forge.jlatexmath.TeXFormula;

/**
 * This class represents the approximator for the grid graphs defined in {@link GridGraph}.
 * <br>
 * The approximator is modelled as a tree.
 * @author Jonas Schulz
 */
public class GridApproximatorTree {
    /**
     * Reference to the root node.
     */
    Node root;
    /**
     * The {@link GridGraph} related to this approximator.
     */
    GridGraph g;
    /**
     * Number of edges in the tree graph visualization.
     */
    int m;
    
    /**
     * An upper limit of the relative estimation error, as maximum ratio of the optimal congestion of b to &Vert; Rb &Vert;<sub>&infin;</sub> for all demands b.
     * This is not well-known, but estimates are discussed in the paper of this thesis.
     * There was no known case that would contradict &alpha; = 3, whereas extensive sampling for d = 1 confirmed that in general,
     * &alpha; has to be at least 3. Still, for the average of those samplings and for all experiments, &alpha; = 2 was a better choice
     * and did not cause any trouble.<br>
     * This parameter can be adjusted with further empiric data.
     */
    double alpha = 3.;
    
    /**
     * Standard constructor.
     * @param g The {@link GridGraph} related to this tree.
     */
    public GridApproximatorTree(GridGraph g){
        int[] lower = new int[g.getDim()], higher = new int[g.getDim()];
        for(int i = 0; i < g.getDim(); i++){
            lower[i] = 0;
            higher[i] = g.nodesPerDim[i]-1;
        }
        root = new Node(g, lower, higher, this);
        this.m = root.m;
        this.g = g;
    }
    
    /**
     * Gets the value of alpha.
     * In case alpha will be set as a formula, this method can be adjusted.
     * @return alpha.
     */
    public double getAlpha(){
//        return Math.log(g.getN())/Math.log(2);
//        return 50.;
        return alpha;
//        return Math.log(g.getN());
    }
    
    /**
     * Updates b for the evaluations related to Rb.
     * @param b Demand to be routed into the approximator.
     */
    public void updateExcessFlows(GridDemand b){
        root.updateExcessFlows(b);
    }
    
    /**
     * Calculates <i>lmax(2&alpha;Rb)</i>.
     * This version does not use shifts.
     * @return <i>lmax(2&alpha;Rb)</i>.
     */
    public double lmax_2alpha_congestion(){
        return Math.log(root.lmax_exp_2alpha_congestion());
    }
    
    /**
     * Calculates e<sup><i>lmax</i>(2&alpha;Rb)</sup> = &sum;<sub>i</sub> e<sup>(2&alpha;Rb)<sub>i</sub></sup> + e<sup>-(2&alpha;Rb)<sub>i</sub></sup>.
     * As in {@link GridFlow}, consider using a shifted version also here to avoid numeric overflow.
     * @return &sum;<sub>i</sub> e<sup>(2&alpha;Rb)<sub>i</sub></sup> + e<sup>-(2&alpha;Rb)<sub>i</sub></sup>.
     */
    public double lmax_exp_2alpha_congestion(){
        return root.lmax_exp_2alpha_congestion();
    }
    
    /**
     * Calculates <i>lmax(2&alpha;Rb)</i>.
     * This version shifts the exponents with the default shift of -2&alpha;&#8729; &Vert;Rb&Vert;<sub>&infin;</sub> 
     * (this way, the maximal exponent after shifting will be 0).
     * The mathematical calculation yields the same as the variant without shifts,
     * but this version prevents overflows.
     * @return <i>lmax(2&alpha;Rb)</i>.
     */
    public double lmax_shifted_2alpha_congestion(){
        double shift = -2*getAlpha()*linf_congestion();
        // alternative:
//        double shift = 2*getAlpha()*linf_congestion() + Math.log(m) - 500;
        return lmax_shifted_2alpha_congestion(shift);
    }
    
    /**
     * Calculates <i>lmax(2&alpha;Rb)</i>.
     * This version shifts the exponents according to the parameter.
     * The mathematical calculation yields the same as the variant without shifts,
     * but this version prevents overflows.
     * @param shift The shift to be <b>added</b> to all exponents.
     * @return <i>lmax(2&alpha;Rb)</i>.
     */
    public double lmax_shifted_2alpha_congestion(double shift){
        return Math.log(lmax_exp_shifted_2alpha_congestion(shift)) - shift;
    }
    
    /**
     * Calculates e<sup><i>lmax</i>(2&alpha;Rb) + <code>shift</code></sup> = &sum;<sub>i</sub> e<sup>(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup> + e<sup>-(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup>.
     * This version uses shifted exponents, according to the parameter <code>shift</code>.
     * @param shift The shift to be <b>added</b> to all exponents.
     * @return &sum;<sub>i</sub> e<sup>(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup> + e<sup>-(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup>.
     */
    public double lmax_exp_shifted_2alpha_congestion(double shift){
        return root.lmax_exp_shifted_2alpha_congestion(shift);
    }
    
    /**
     * This inner class represents a node of the approximator tree structure.<br>
     * Each node can be interpreted as a hypercube cut of the original graph.
     * Given a demand <i>b</i>, the optimal congestion of a flow that satisfies <i>b</i>
     * can not be better than b<sub>C</sub>/c<sub>C</sub>, where C is a cut, c<sub>C</sub> its capacity and
     * b<sub>C</sub> the sum of all <i>b<sub>i</sub></i> for the vertices <i>i</i> in C.
     * This is the metric represented in each (Rb)<sub>i</sub>, and each node calculates one of those metrics.
     */
    static class Node {
        /**
         * Lower bounds of the hypercube.
         */
        int[] lowerIndices;
        /**
         * Higher bounds of the hypercube.
         */
        int[] higherIndices;
        /**
         * Capacity of the hypercube cut.
         * The hypercube is defined by {@link #lowerIndices} and {@link #higherIndices}.
         */
        int capacity_cut;
        /**
         * Current excess flow of this hypercube.
         * The hypercube is defined by {@link #lowerIndices} and {@link #higherIndices}.
         * The excess flow has to be updated via {@link #updateExcessFlows(GridDemand b)}.
         */
        double current_excess_flow;
        /**
         * The children of this {@link Node} represent the hypercube splits from {@link GridGraph#split(int[], int[])}.
         * If this {@link Node} represents a single hypercube, children is {@code null}.
         */
        Node[] children;
        /**
         * Reference to the parent node. Currently not in use.
         */
        Node parent;
        /**
         * Current value to be used for <i>&nabla;lmax(2&alpha;Rb)</i>.
         * The gradient has to be updated via {@link #set_edge_gradient_2alpha_potential()}, {@link #set_edge_gradient_shift_2alpha_potential()} or
         * {@link #set_edge_gradient_shift_2alpha_potential(double)}. Doing so will set the gradient according to <code>b</code> from the last call of 
         * {@link #updateExcessFlows(GridDemand)}.
         */
        double current_edge_grad;
        /**
         * The {@link GridApproximatorTree} object maintaining the tree structure with this {@link Node}.
         */
        GridApproximatorTree T;
        /**
         * The number of edges in the subtree where this node is the root.
         */
        int m;

        /**
         * Recursive constructor.
         * @param g The grid graph on which the approximator operates.
         * @param lowerIndices Lower bounds of the hypercube represented by this node.
         * @param higherIndices Higher bounds of the hypercube represented by this node.
         * @param tree {@link GridApproximatorTree} object of the tree in which this node lies.
         */
        public Node(GridGraph g, int[] lowerIndices, int[] higherIndices, GridApproximatorTree tree){
            this.lowerIndices = lowerIndices;
            this.higherIndices = higherIndices;
            this.capacity_cut = g.capHyperBox(lowerIndices, higherIndices);
            int dif = g.different(lowerIndices, higherIndices);
            this.T = tree;
            this.m = 0;
            if(dif == 0){
                this.children = null;
            } else {
                int[][][] splits = g.split(lowerIndices, higherIndices);
                this.children = new Node[splits.length];
                for(int i = 0; i < splits.length; i++){
                    children[i] = new Node(g, splits[i][0], splits[i][1], T);
                    children[i].parent = Node.this;
                    this.m += 1 + children[i].m;
                }
            }
            this.current_excess_flow = 0;
        }

        /**
         * Standard setter for {@link #current_excess_flow}.
         * @param current_excess_flow New value for {@link #current_excess_flow}.
         */
        public void setCurrent_excess_flow(double current_excess_flow) {
            this.current_excess_flow = current_excess_flow;
        }
        
        /**
         * Returns whether this node is a leaf. A leaf represents only a single vertex.
         * @return <code>true</code> iff this node is a leaf.
         */
        public boolean isLeaf(){
            return children == null;
        }
        
        /**
         * Returns whether this node is the root of {@link #T}.
         * The root represents the whole graph, and thus has no entry in Rb.
         * @return <code>true</code> iff this node is the root of {@link #T}.
         */
        public boolean isRoot(){
            return this.capacity_cut == 0;
        }
        
        /**
         * Sets the leaves to the given excess flows b and recursively (post-order) calculates the
         * excess flows of higher-level hypercube cuts of the graph, that is, the inner
         * nodes of this tree structure.
         * The congestion approximation at a node v, (Rb)<sub>v</sub>, directly can be read via 
         * <code>(v.current_excess_flow/v.capacity_cut)</code>.
         * @param b Demand to be routed into R. Note that you often want to route some residual demand here, i.e. <i>(b-Bf)</i>.
         * @return Excess flow of this subtree, for recursive use.
         */
        public double updateExcessFlows(GridDemand b){
            if(isLeaf()){
                // int index = T.g.toIndex(lowerIndices);
                this.current_excess_flow = b.get(lowerIndices);
            } else {
                double sum = 0;
                for(Node n : children){
                    sum += n.updateExcessFlows(b);
                }
                this.current_excess_flow = sum;
            }
            return this.current_excess_flow;
        }
        
        /**
         * Calculates e<sup>lmax(2&alpha;Rb)</sup> for the parts of R and b in this sub-tree.
         * The concrete calculation is &sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub></sup>).
         * Set b via updateExcessFlows(b).
         * The exponent level is kept to optimize calculation efficiency.
         * To get lmax(2&alpha;Rb), simply pipe this result in Math.log or
         * use the {@link #lmax_2alpha_congestion()}, {@link #lmax_shifted_2alpha_congestion()} or {@link #lmax_shifted_2alpha_congestion(double)}
         * methods of {@link #T}.
         * @return &sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub></sup>).
         */
        public double lmax_exp_2alpha_congestion(){
            double sum = 0;
            if(!isRoot()){
                double cong = this.current_excess_flow / this.capacity_cut;
                cong *= 2*T.getAlpha();
                double exp = Math.exp(cong);
                sum += exp + (1./exp);
            }
            if(!isLeaf()){
                for(Node n : children){
                    sum += n.lmax_exp_2alpha_congestion();
                }
            }
            return sum;
        }
        
        /**
         * Similar to {@link #lmax_exp_2alpha_congestion()}, but with a shift to avoid numeric instability.
         * Actual calculation is &sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub> + <code>shift</code></sup>+exp<sup>-(2&alpha;Rb)<sub>i</sub> + <code>shift</code></sup>).<br>
         * This yields e<sup>lmax(2&alpha;Rb)+<code>shift</code></sup>.
         * @param shift Shift to be <b>added</b> to the exponents.
         * @return &sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub> + <code>shift</code></sup>+exp<sup>-(2&alpha;Rb)<sub>i</sub> + <code>shift</code></sup>).
         */
        public double lmax_exp_shifted_2alpha_congestion(double shift){
            double sum = 0;
            if(!isRoot()){
                double cong = this.current_excess_flow / this.capacity_cut;
                cong *= 2*T.getAlpha();
                double exp = Math.exp(cong + shift);
                sum += exp + Math.exp((-cong) + shift);
            }
            if(!isLeaf()){
                for(Node n : children){
                    sum += n.lmax_exp_shifted_2alpha_congestion(shift);
                }
            }
            return sum;
        }

        /**
         * Calculates <i>&sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub></sup>)&#8729;&nabla;lmax(2&alpha;Rb)<sub>j</sub></i> for the respective j of this node.
         * The concrete calculation is (e<sup>(2&alpha;Rb)<sub>j</sub></sup>-e<sup>-(2&alpha;Rb)<sub>j</sub></sup>).<br>
         * Set b via {@link #updateExcessFlows(GridDemand b)}.
         * The term <i>&sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub></sup>)</i> evaluates to the same
         * value for all j. Hence, some linear maps on <i>&nabla;lmax</i> can be accelerated by postponing this factorization.<br>
         * The actual gradient is currently never used directly; this method implicitly stores it inside this object, whereas other
         * functions involving the gradient (e.g. <i>l<sub>&infin;</sub></i>) are implemented in this class and return the respective concrete mathematical objects.
         */
        public void set_edge_gradient_2alpha_potential_rec(){
            // only calculating the non-constant term e^xj - e^(-xj)
            if(!isRoot()){
                double exp = Math.exp(2*T.getAlpha()*this.current_excess_flow / this.capacity_cut);
                this.current_edge_grad = exp - (1./exp);
            }
            if(!isLeaf()){
                for(Node n : children){
                    n.set_edge_gradient_2alpha_potential_rec();
                }
            }
        }
        /**
         * Same as {@link #set_edge_gradient_2alpha_potential_rec()}, but using shifts for the exponents
         * to prevent numeric overflows.
         * @param shift The shift to be <b>added</b> to all exponents.
         */
        public void set_edge_gradient_shift_2alpha_potential_rec(double shift){
            // only calculating the non-constant term e^xj - e^(-xj)
            if(!isRoot()){
                double a2c = 2*T.getAlpha()*this.current_excess_flow / this.capacity_cut;
                double exp = Math.exp(a2c + shift);
                double exp2 = Math.exp((-a2c) + shift);
                this.current_edge_grad = exp - exp2;
            }
            if(!isLeaf()){
                for(Node n : children){
                    n.set_edge_gradient_shift_2alpha_potential_rec(shift);
                }
            }
        }
        
        /**
         * Calculates this node's contribution to <i>R<sup>T</sup>&#8729; &nabla;lmax(2&alpha;Rb)</i>,
         * as it can be expressed as <i>&sum;<sub>i</sub> (g<sub>i</sub> / c<sub>i</sub>)</i> &#8729; I<sub>i</sub>, where
         * <ul>
         * <li><i>g<sub>i</sub></i> is the gradient for node i,</li>
         * <li><i>c<sub>i</sub></i> is the capacity of the cut represented by node i and</li>
         * <li><i>I<sub>i</sub></i> is the incidence vector for node <i>i</i> (i.e., (1,0,0,1) iff graph nodes 1 and 4 are in this box).</li>
         * </ul>
         * This method recursively calls all children, s.t. calling the root node will store <i>R<sup>T</sup>&#8729; &nabla;lmax(2&alpha;Rb)</i> in <code>d</code>.
         * @param d {@link GridDemand} to store the results in.
         */
        public void add_edge_gradients_to_incident_node_excesses_v2(GridDemand d){
            if(!isRoot()){
                int[] indices = T.g.indicesOfBoxNodes(lowerIndices, higherIndices);
                for(int i : indices){
                    d.add(i, current_edge_grad/this.capacity_cut);
                }
            }
            if(!isLeaf()){
                for(Node n : children){
                    n.add_edge_gradients_to_incident_node_excesses_v2(d);
                }
            }
        }
        
        /**
         * {@inheritDoc } 
         */
        @Override
        public String toString() {
            return leadingWhitespaceString(0);
        }
        
        /**
         * Implements functionality of the {@link #toString()} method.
         * Recursively constructs the string representing this whole subtree.
         * @param whitespaces Extra whitespace indent per subtree depth.
         * @return String representation of this subtree to be embedded into the whole tree string representation.
         */
        private String leadingWhitespaceString(int whitespaces){
            String lws = "";
            for(int i = 0; i < whitespaces; i++){ lws += " "; }
            String s = lws + (isLeaf() ? "Leaf:\n" : "Node:\n");
            s += lws + "Lower Indices: "+Arrays.toString(lowerIndices)+"\n";
            s += lws + "Higher Indices: "+Arrays.toString(higherIndices)+"\n";
            s += lws + "Capacity: "+capacity_cut+"\n";
            s += lws + "Current Excess Flow: "+current_excess_flow+"\n";
            s += lws + "Current Gradient: "+this.current_edge_grad+"\n";
            if(!isLeaf()){
                for(Node n : children){
                    s += n.leadingWhitespaceString(whitespaces + 4);
                }
            }
            return s;
        }

        /**
         * Calculates the <i>l<sub>&infin;</sub></i>-norm of Rb.
         * @return &Vert;Rb&Vert;<sub>&infin;</sub>
         */
        double linf_congestion() {
            double max = 0.;
            if(!isRoot()){
                max = Math.abs(current_excess_flow / capacity_cut);
            }
            if(!isLeaf()){
                for(Node n : children){
                    double d = n.linf_congestion();
                    if(d > max)
                        max = d;
                }
            }
            return max;
        }

        /**
         * Returns the height of this subtree.
         * @return The height of this subtree. Leaves have height 0.
         */
        int height(){
            if(isLeaf()) return 0;
            int max = 0;
            for(Node n : children){
                if(n.height() > max) max = n.height();
            }
            return 1 + max;
        }
        
        /**
        * Returns TikZ-input for a visual representation of this subtree, including the current values of <i>lmax(2&alpha;Rb)</i> and <i>&nabla;lmax(2&alpha;Rb)</i>.<br>
        * @param xpos x-position of this node. Recursive calls adjust this parameter with the volume of the respective hypercube represented by each child.
        * @return TikZ input for an extended visual representation of this subtree.
        */
        private String tikz2D(double xpos) {
            int yscale = 4;
            int ypos = yscale*height();
            DecimalFormat df = new DecimalFormat("#.0000");
            DecimalFormat gf = new DecimalFormat("#0.##E0");
//            int xpos = T.g.volume(lowerIndices, higherIndices)/2;
//            xpos += T.g.toIndex(lowerIndices);
            String s = "";
            String ns = isLeaf() ? ""+T.g.toIndex(lowerIndices) : "";
            s += "\\node[roundnode, minimum size = 0.9cm] at ("+(xpos+T.g.volume(lowerIndices, higherIndices)/2.)+", "+ypos+") ("+((int)xpos)+"h"+ypos+") {"+ns+"};\n";
            if(!isLeaf()){
                double sum = xpos;
                for(Node n : children){
                    String col = (n.current_excess_flow == 0.) ? "black!20" : "black";
                    String colGrad = (n.current_edge_grad == 0.) ? "red!20" : "red";
                    String colEdge = (n.current_excess_flow == 0.) ? "black!20" : "blue!40";
                    s += n.tikz2D(sum + T.g.volume(n.lowerIndices, n.higherIndices)/2.);
                    s += String.format("\\draw[thick, ->, color=%s] (",colEdge)+(int)xpos+"h"+ypos+") -- ("+((int)(sum + T.g.volume(n.lowerIndices, n.higherIndices)/2.))+"h"+(yscale*n.height())+")"
                            + String.format(" node[sloped,midway,above=-0.1cm] {\\textcolor{%s}{%s/%d}, \\textcolor{%s}{%s}}", col, df.format(n.current_excess_flow), n.capacity_cut, colGrad, gf.format(n.current_edge_grad)) + ";\n";
                    sum += T.g.volume(n.lowerIndices, n.higherIndices);
                }
            }
            return s;
        }
    }

    /**
     * {@inheritDoc }
     */
    @Override
    public String toString() {
        String s = "";
        s += ("Grid Dimensions: "+Arrays.toString(g.nodesPerDim)) + "\n";
        s += ("Root:") + "\n";
        s += (root) + "\n";
        return s;
    }
    
    /** 
     * Calculate <i>&sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub></sup>)&#8729;&nabla;lmax(2&alpha;Rb)</i> as gradient of tree edges,
     * see the calculation in [She13] of <i>p<sub>2</sub> = &nabla;lmax(x<sub>2</sub>)</i> with <i>x<sub>2</sub> = 2&alpha; &#8729; R (b - Bf)</i>.
     * For more implementation details, see the documentation of {@link Node#set_edge_gradient_2alpha_potential_rec() }.
     * The <i>1/(&sum;<sub>i</sub>e<sup>(2&alpha;Rb)<sub>i</sub></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub></sup>)</i> part is omitted, as there are only linear 
     * operations, thus this scalar can be factored out, postponing the multiplication and evaluation.
     */
    public void set_edge_gradient_2alpha_potential(){
        // if not factoring out, add following line:
        // double expsum = root.lmax_exp_2alpha_congestion();
        root.set_edge_gradient_2alpha_potential_rec();
    }
    
    /**
     * Similar to {@link #set_edge_gradient_2alpha_potential()}, but shifts the exponents.
     * The concrete formula is then <i>&sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup>)&#8729;&nabla;lmax<sub><code>shift</code></sub>(2&alpha;Rb)</i>
     * with <i>(&nabla;lmax<sub><code>shift</code></sub>(2&alpha;Rb))<sub>j</sub> = (e<sup>(2&alpha;Rb)<sub>j</sub>+<code>shift</code></sup> - e<sup>-(2&alpha;Rb)<sub>j</sub>+<code>shift</code></sup>)
     *  / &sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup>)</i>.
     * @param shift 
     */
    public void set_edge_gradient_shift_2alpha_potential(double shift){
        root.set_edge_gradient_shift_2alpha_potential_rec(shift);
    }
    
    /**
     * Same as {@link #set_edge_gradient_shift_2alpha_potential(double)}, but with automatic shift s.t. all exponents will be at most 0.
     */
    public void set_edge_gradient_shift_2alpha_potential(){
        double shift = getDefaultShift();
        set_edge_gradient_shift_2alpha_potential(shift);
    }
    
    /**
     * Calculates <i>s<sub>b</sub> &#8729; R<sup>T</sup> &#8729; &nabla;lmax(2&alpha;Rb)</i> with 
     * <i>s<sub>b</sub> = &sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub></sup>)&#8729;&nabla;lmax(2&alpha;Rb)</i>.
     * Use {@link #set_edge_gradient_2alpha_potential()}, {@link #set_edge_gradient_shift_2alpha_potential()} or {@link #set_edge_gradient_shift_2alpha_potential(double)}
     * to set <i>b</i> first.<br>
     * If the shifted variants were used, the returned demand is instead <i>sh<sub>b</sub> &#8729; R<sup>T</sup> &#8729; &nabla;lmax<sub><code>shift</code></sub>(2&alpha;Rb)</i>
     * with <i>sh<sub>b</sub> = &sum;<sub>i</sub>(e<sup>(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup>+e<sup>-(2&alpha;Rb)<sub>i</sub>+<code>shift</code></sup>)</i>.
     * @return <i>s<sub>b</sub> &#8729; R<sup>T</sup> &#8729; &nabla;lmax(2&alpha;Rb)</i>, or <i>sh<sub>b</sub> &#8729; R<sup>T</sup> &#8729; &nabla;lmax<sub><code>shift</code></sub>(2&alpha;Rb)</i>
     * when shifts were used at gradient calculation.
     */
    public GridDemand mult_Rt_edge_gradient(){
        GridDemand d = new GridDemand(this.g);
        root.add_edge_gradients_to_incident_node_excesses_v2(d);
        return d;
    }
    
    /**
     * Calculates the default shift.
     * The default shift shifts all exponents s.t. the maximum exponent will be zero.
     * As the exponents are <i>&plusmn;(2&alpha;Rb)<sub>i</sub></i>, the shift is thus chosen as
     * <i>-2&alpha; &#8729; &Vert;Rb&Vert;<sub>&infin;</sub></i>.
     * @return The default shift.
     */
    public double getDefaultShift(){
        return -2*getAlpha()*linf_congestion();
    }
    
    /**
     * Calculates some useful metric of [She13].
     * Make sure that the current residual demand, <i>r<sup>(i)</sup> = (b-Bf<sup>(i)</sup>)</i> is routed in the gradient, not b!
     * @param b Original demand vector.
     * @return <i>b<sup>T</sup>v / &Vert;B<sup>T</sup>v&Vert;<sub>1</sub></i> with <i>v</i> as result of {@link #mult_Rt_edge_gradient()}.
     */
    public double compute_lemma_congestion(GridDemand b){
        GridDemand v = mult_Rt_edge_gradient();
        double z = GridDemand.scalar_prod(b, v);
        double n = v.toPotentialDiffEdgesFlow().l1();
        return z/n;
    }
    
    /**
     * Calculates some useful metric of [She13].
     * It is calculated as <i>(1+&epsilon;) &#8729; {@link #compute_lemma_congestion(GridDemand b)}</i>.
     * @param lemma_congestion Use this to input the value you got from {@link #compute_lemma_congestion(GridDemand b)}.
     * @param eps The value of &epsilon;.
     * @return The upper limit from the referred metric in [She13].
     */
    public double compute_lemma_congestion_upper(double lemma_congestion, double eps){
        return (1+eps)*lemma_congestion;
    }
    
    /**
     * Calculates <i>&Vert;Rb&Vert;<sub>&infin;</sub></i> for the current excess flow b.
     * Use {@link #updateExcessFlows(GridDemand b)} to set b.
     * @return <i>&Vert;Rb&Vert;<sub>&infin;</sub></i>.
     */
    public double linf_congestion(){
        return root.linf_congestion();
    }
    
    /**
     * Returns TikZ-input for a visual representation of this approximator, including the current values of <i>lmax(2&alpha;Rb)</i> and <i>&nabla;lmax(2&alpha;Rb)</i>.<br>
     * @return TikZ input for an extended visual representation of this approximator.
     */
    public String tikz2D(){
        String prefix = "\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]\n";
        prefix += "\\hspace{-4cm}\n";
        String suffix = "\\end{tikzpicture}";
        return prefix+root.tikz2D(0)+suffix;
    }
    
    /**
     * Renders the TikZ-input of {@link #tikz2D()}  into an image.
     * @deprecated Not implemented.
     */
    @Deprecated
    public void createPNG(){
//        TeXFormula lat = new TeXFormula(); // this line alone leads to errors.
//        getClass().getClassLoader().getResourceAsStream("PNGs/");
//        getClass().getClassLoader().getResourceAsStream("PNGs\\");
//        String firstLine = "\\hspace{-4cm}";
//        firstLine += String.format("\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=%1.1f]", 1.2*(g.getN()));
//        String lastLine = "\\end{tikzpicture}";
//        TeXFormula formula = new TeXFormula(firstLine + tikz2D() + lastLine);
//        formula.createPNG(TeXConstants.STYLE_DISPLAY, 100, "PNGs/Tree-"+System.currentTimeMillis(), Color.white, Color.black);
    }
    
}
