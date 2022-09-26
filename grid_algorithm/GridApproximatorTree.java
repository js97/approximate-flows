package grid_algorithm;

import java.awt.Color;
import java.text.DecimalFormat;
import java.util.Arrays;

import org.scilab.forge.jlatexmath.TeXConstants;
import org.scilab.forge.jlatexmath.TeXFormula;

/**
 *
 * @author kroka
 */
public class GridApproximatorTree {
    Node root;
    GridGraph g;
    int m;
    
    double alpha = 10.;
    
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
    public double getAlpha(){
//        return Math.log(g.getN())/Math.log(2);
//        return 50.;
        return alpha;
//        return Math.log(g.getN());
    }
    
    public void updateExcessFlows(GridDemand b){
        root.updateExcessFlows(b);
    }
    public double lmax_2alpha_congestion(){
        return Math.log(root.lmax_exp_2alpha_congestion());
    }
    public double lmax_exp_2alpha_congestion(){
        return root.lmax_exp_2alpha_congestion();
    }
    public double lmax_shifted_2alpha_congestion(){
        double shift = -2*getAlpha()*linf_congestion();
        // alternative:
//        double shift = 2*getAlpha()*linf_congestion() + Math.log(m) - 500;
        return lmax_shifted_2alpha_congestion(shift);
    }
    public double lmax_shifted_2alpha_congestion(double shift){
        return Math.log(lmax_exp_shifted_2alpha_congestion(shift)) - shift;
    }
    public double lmax_exp_shifted_2alpha_congestion(double shift){
        return root.lmax_exp_shifted_2alpha_congestion(shift);
    }
    
    static class Node {
        int[] lowerIndices;
        int[] higherIndices;
        int capacity_cut;
        double current_excess_flow;
        Node[] children;
        Node parent;
        double current_edge_grad;
        //TODO implement parent functionality
        GridApproximatorTree T;
        int m;

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
//        public Node(int[] lowerIndices, int[] higherIndices, int capacity_cut) {
//            this.lowerIndices = lowerIndices;
//            this.higherIndices = higherIndices;
//            this.capacity_cut = capacity_cut;
//            this.current_excess_flow = 0;
//            this.children = null;
//        }

        public void setCurrent_excess_flow(double current_excess_flow) {
            this.current_excess_flow = current_excess_flow;
        }
        public boolean isLeaf(){
            return children == null;
        }
        public boolean isRoot(){
            return this.capacity_cut == 0;
        }
        
        /**
         * Sets the leaves to the given excess flows b and calculates the
         * excess flows of higher-level parts of the graph, that is, the inner
         * nodes of this tree structure.
         * The congestion approximation at a node v, (Rb)_v, directly can be read via 
         * v.current_excess_flow/v.capacity_cut.
         * @param b
         * @return 
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
         * Calculates exp(lmax(2 alpha * R b)) for the parts of R and b in this sub-tree.
         * The concrete calculation is sum(exp((2 alpha * Rb)_i)+exp(-(2 alpha * Rb)_i)).
         * Set b via updateExcessFlows(b).
         * The exponent level is kept to optimize calculation efficiency.
         * To get lmax(2 alpha * R b), simply pipe this result in Math.log or
         * use the lmax_2alpha_congestion() or lmax_shifted_2alpha_congestion()
         * methods of the whole tree data structure.
         * @return 
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
         * Same as lmax_exp_2alpha_congestion, but with a shift to avoid numeric instability.
         * Actual calculation is sum(exp((2 alpha * Rb)_i + shift)+exp(-(2 alpha * Rb)_i) + shift).
         * @param shift
         * @return 
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
         * Calculates this node's contribution to Rt*grad(lmax(2 alpha R b_residual)),
         * as it can be expressed as sum(i over nodes) (g_i * I_i), where
         * - g_i is the gradient for node i
         * - I_i is the incidence vector for node i (i.e., (1,0,0,1) iff graph nodes 1 and 4 are in this box)
         * @param d demand vector to store results in
         * @deprecated Capacity bug. Use fixed version v2 instead.
         */
        @Deprecated
        public void add_edge_gradients_to_incident_node_excesses(GridDemand d){
            if(!isRoot()){
                int[] indices = T.g.indicesOfBoxNodes(lowerIndices, higherIndices);
                for(int i : indices){
                    d.add(i, current_edge_grad);
                }
            }
            if(!isLeaf()){
                for(Node n : children){
                    n.add_edge_gradients_to_incident_node_excesses(d);
                }
            }
        }
        public void add_edge_gradients_to_incident_node_excesses_v2(GridDemand d){
            if(!isRoot()){
                int[] indices = T.g.indicesOfBoxNodes(lowerIndices, higherIndices);
//                System.out.println(String.format("Adding %f/%d to indices %s", current_edge_grad, capacity_cut, Arrays.toString(indices)));
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
        
        @Override
        public String toString() {
            return leadingWhitespaceString(0);
        }
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
        @Deprecated
        double org_demand;
        @Deprecated
        /**
         * Not implemented yet.
         */
        double demand_bound() {
            double pot = 0.;
            double orgdem = 0.;
            if(isLeaf()) return pot * orgdem;
            double sum = 0;
            for(Node n : children){
                sum += n.demand_bound();
            }
            return sum;
        }
        @Deprecated
        /**
         * Debug Purpose.
         * @param b 
         */
        void set_orgdemand(GridDemand b){
            if(isLeaf()){
                // int index = T.g.toIndex(lowerIndices);
                this.org_demand = b.get(lowerIndices);
                return;
            } else {
                for(Node n : children){
                    n.set_orgdemand(b);
                }
            }
        }

        int height(){
            if(isLeaf()) return 0;
            int max = 0;
            for(Node n : children){
                if(n.height() > max) max = n.height();
            }
            return 1 + max;
        }
        
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

    @Override
    public String toString() {
        String s = "";
        s += ("Grid Dimensions: "+Arrays.toString(g.nodesPerDim)) + "\n";
        s += ("Root:") + "\n";
        s += (root) + "\n";
        return s;
    }
    
    /** Calculate grad(2 alpha * R * (excess flows)) as gradient of tree edges,
    * see calculation in paper of p_2 = grad(lmax(x_2)).
    * The 1/sum(exp(xi)+exp(-xi)) part is omitted, as there are only linear 
    *    operations, thus this scalar can be factored out.
    * */
    public void set_edge_gradient_2alpha_potential(){
        // if not factoring out, add following line:
        // double expsum = root.lmax_exp_2alpha_congestion();
        root.set_edge_gradient_2alpha_potential_rec();
    }
    public void set_edge_gradient_shift_2alpha_potential(double shift){
        root.set_edge_gradient_shift_2alpha_potential_rec(shift);
    }
    public void set_edge_gradient_shift_2alpha_potential(){
        double shift = getDefaultShift();
        set_edge_gradient_shift_2alpha_potential(shift);
    }
    public GridDemand mult_Rt_edge_gradient(){
        GridDemand d = new GridDemand(this.g);
//        root.add_edge_gradients_to_incident_node_excesses(d);
        root.add_edge_gradients_to_incident_node_excesses_v2(d);
        return d;
    }
    public double getDefaultShift(){
//        System.out.println("Calculating shift... tree:\n"+this.tikz2D());
//        System.out.println("Shift will be "+(-2*getAlpha()*linf_congestion()));
        return -2*getAlpha()*linf_congestion();
//        return -linf_congestion();
    }
    
    /**
     * Make sure that (b-Bf) is routed in gradient, not b!
     * @param b Since b is not explicitly required for calculating maxflows, but
     * required for calculating this useful metric, it is requested as input.
     * @return 
     */
    public double compute_lemma_congestion(GridDemand b){
        GridDemand v = mult_Rt_edge_gradient();
        double z = GridDemand.scalar_prod(b, v);
        double n = v.toPotentialDiffEdgesFlow().l1();
        return z/n;
    }
    public double compute_lemma_congestion_upper(double lemma_congestion, double eps){
        return (1+eps)*lemma_congestion;
    }
    
    /**
     * Calculates ||Rb||_inf for the current excess flow b.
     * Use updateExcessFlows(b) to set b.
     * @return ||Rb||_inf.
     */
    public double linf_congestion(){
        return root.linf_congestion();
    }
    
    public String tikz2D(){
        String prefix = "\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]\n";
        prefix += "\\hspace{-4cm}\n";
        String suffix = "\\end{tikzpicture}";
        return prefix+root.tikz2D(0)+suffix;
    }
    
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
