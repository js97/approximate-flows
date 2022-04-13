/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
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
return 50.;
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
         */
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
        root.add_edge_gradients_to_incident_node_excesses(d);
        return d;
    }
    public double getDefaultShift(){
        return -2*getAlpha()*linf_congestion();
    }
    
    public double linf_congestion(){
        return root.linf_congestion();
    }
    
    public String tikz2D(){
        return root.tikz2D(0);
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
