/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grid_algorithm;

import com.opencsv.CSVWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map.Entry;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;
import project_utils.DoubleSequence;
import project_utils.Tuple;

/**
 *
 * @author kroka
 */
public class GridApproximation {
    GridGraph g;
    GridApproximatorTree t;
    double currentScale = 1.;
    
    DoubleSequence stepSize = DoubleSequence.one;
    
    public GridApproximation(GridGraph g){
        // debug stuff
        debug_csv : {
            String subfolder = "";
            String path = "C:/Users/kroka/Documents/Masterarbeit/Logs/";
            try {
                this.writer = new CSVWriter(new FileWriter(path + subfolder + "deltas.csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);
                writer.writeNext(new String[]{"delta"});
                this.iterationWriter = new CSVWriter(new FileWriter(path + subfolder + "potentials.csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);
                // <editor-fold defaultstate="collapsed" desc="Some code">
                iterationWriter.writeNext(new String[]{"potential", "current scale", String.format("edgeflow-%d-%d",4,8), "edgeflow-rescaled"});
                // </editor-fold>
            } catch (IOException ex) {
                Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        this.g = g;
        this.t = new GridApproximatorTree(g);
    }
    public GridApproximation(GridGraph g, String input_identifier){
        // debug stuff
        debug_csv : {
            String subfolder = "";
            String path = "C:/Users/kroka/Documents/Masterarbeit/Logs/";
            try {
                this.writer = new CSVWriter(new FileWriter(path + subfolder + input_identifier + "_deltas.csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);
                writer.writeNext(new String[]{"delta"});
                this.iterationWriter = new CSVWriter(new FileWriter(path + subfolder + input_identifier + "_potentials.csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);
                iterationWriter.writeNext(new String[]{"potential", "current scale", String.format("edgeflow-%d-%d",4,8), "edgeflow-rescaled"});
            } catch (IOException ex) {
                Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        this.g = g;
        this.t = new GridApproximatorTree(g);
    }
    
    boolean printTikz = true;
    boolean printToString = false;
    int iterationLimit = 515;
    // debug purpose
    CSVWriter writer, iterationWriter;
    
    public Tuple<HashedGridFlow, HashedGridFlow> AlmostRoute_inputs(GridDemand b, double eps, double alpha, DoubleSequence h){
        t.alpha = alpha;
        stepSize = h;
        return AlmostRoute(b, eps);
    }
    
    // just to record how many iterations have been made
    int iterations = 0;
    public Tuple<HashedGridFlow, HashedGridFlow> AlmostRoute(GridDemand b, double eps){
        // initialization
        t.updateExcessFlows(b);
        double linf = t.linf_congestion();
        double s = ((16/eps)*Math.log(g.getN()))/(2*t.getAlpha()*linf);
        b = b.scale(s);
        currentScale = s;
        HashedGridFlow currentFlow = new HashedGridFlow(g);
        Tuple<HashedGridFlow, HashedGridFlow> iter_result;
        double delta;
        // iteration
        do{
            // [bullet 1] and pre-calculation phi, gradphi
//            System.out.println("Current Tree (before iter):");
//            printTree(t);
            iter_result = iteration(currentFlow, b, eps);
            // [bullet 2] delta calculation
            delta = iter_result.b.l1();
//            System.out.println("Current scale: "+this.currentScale);
            // [bullet 3] update flow approximation
            if(delta >= eps/4){
                double val = (-1)*(delta/(1+4*t.getAlpha()*t.getAlpha()));
                // scale step with stepSize (experimental):
                val *= stepSize.at(iterations);
                for(Entry<Integer, Double> e : iter_result.b.entries.entrySet()){
                    // capacity is 1 for all edges here
//                    double val = -(delta/(1+4*t.getAlpha()*t.getAlpha()))*e.getValue();
                    currentFlow.add(e.getKey(), val * Math.signum(e.getValue()));
                }
            }
            // just for information
            iterations++;
//            System.out.println(String.format("Delta = %f", delta));
            debug_csv : {
                writer.writeNext(new String[]{""+delta});
            }
            // loop = "repeat" (halting condition = [bullet 4])
        } while (delta >= eps/4 && iterations < iterationLimit);
        // return
        iter_result.a = iter_result.a.scale(1./currentScale);
        iter_result.b = iter_result.b.scale(1./currentScale);
//debug : {
//    System.out.println("Scaled Demand:\\\\");
//    boolean tmp = printTikz;
//    printTikz = true;
//    printDemand(b);
//    double sdb = 1./currentScale;
//    b = b.scale(sdb);
//    System.out.println(String.format("Re-scaled Demand with total scale %f:\\\\", currentScale));
//    printDemand(b);
//    printTikz = tmp;
//}
        System.out.println("Total scaling: "+currentScale);
        System.out.println("Stopped at iteration "+iterations);
        
        debug_csv : {
            try {
                writer.flush();
                iterationWriter.flush();
            } catch (IOException ex) {
                Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        return iter_result;
    }
    
    public HashedGridFlow CompleteRoute_inputs(GridDemand b, double eps, double alpha, DoubleSequence h){
        t.alpha = alpha;
        stepSize = h;
        return CompleteRoute(b, eps);
    }
    public HashedGridFlow CompleteRoute(GridDemand b, double eps){
        // first route: AlmostRoute
        Tuple<HashedGridFlow, HashedGridFlow> f0_ = AlmostRoute(b, eps);
        int m = g.getM();
        int T = (int)Math.log(2*m);
        GridDemand bi = b;
        HashedGridFlow fi = f0_.a;
        for(int i = 1; i <= T; i++){
            // consecutive routings: AlmostRoute with residual and eps=0.5
            bi = GridDemand.subtract(bi, fi.calculateExcessFlows());
            Tuple<HashedGridFlow, HashedGridFlow> fi_ = AlmostRoute(bi, 0.5);
            fi = fi_.a;
        }
        //TODO: Spanning Tree routing
        return fi;
    }
    
    
    // part of the paper's "repeat" part, namely the scaling of f and b (bullet 1),
    // and the calculation of potential and potential gradient of f (to be used for other bullet points)
    private Tuple<HashedGridFlow, HashedGridFlow> iteration(HashedGridFlow currentFlow, GridDemand demand, double eps){
        double pot = potential(currentFlow, demand);
//        System.out.println("Pot @ "+iterations+" :   "+pot);
        double locScale = 1.;
//        System.out.println(String.format("Comparing %f with %f...",pot,16*(1./eps)*Math.log(g.getN())));
        while(pot < 16*(1./eps)*Math.log(g.getN())){
            // inplace scaling since caller method relies on referred objects
            currentFlow = currentFlow.scale_inplace(17./16);
            demand.scale_inplace(17./16);
            this.currentScale *= 17./16;
            locScale *= 17./16;
            pot = potential(currentFlow, demand);
        }
        debug_csv : {
            double edgeFlow = currentFlow.get(4, 8);
            double edgeFlowRescaled = edgeFlow / currentScale;
            iterationWriter.writeNext(new String[]{""+pot, ""+currentScale, ""+edgeFlow, ""+edgeFlowRescaled});
        }
//        System.out.println("Scaled Potential: "+pot);
//        System.out.println("         > Graph: "+cpg);
//        System.out.println("         >  Tree: "+cpt);
//        System.out.println("     >  Local Scale: "+locScale);
//        System.out.println("     >  Total Scale: "+currentScale);
        // GridFlow grad_pot = new GridFlow(g);
        HashedGridFlow grad_potential = grad_potential(currentFlow, demand);
//        System.out.println("Current Flow: ");
//        printFlow(currentFlow);
//        System.out.println("Current Gradient: ");
//        printFlow(grad_potential);
        
        
        return new Tuple<>(currentFlow, grad_potential);
    }
    double cpt = 0., cpg = 0.;
    public double potential(HashedGridFlow currentFlow, GridDemand demand){
//        System.out.println(String.format("Calculating Potential in iteration %d....",this.iterations));
//        System.out.println("Given Flow:");
//        printFlow(currentFlow);
//        System.out.println(currentFlow);
//        System.out.println("Given Demand: ");
//        printDemand(demand);
        GridDemand bf = currentFlow.calculateExcessFlows();
//        System.out.println("Calculated Excess Flows:");
//        printDemand(bf);
        GridDemand residualDemand = GridDemand.subtract(demand, bf);
//        System.out.println("Calculated Residual Demand: ");
//        printDemand(residualDemand);
        double pot_graph = currentFlow.lmax_shifted();
        t.updateExcessFlows(residualDemand);
        double pot_tree = t.lmax_shifted_2alpha_congestion();
        this.cpt = pot_tree;
        this.cpg = pot_graph;
//        System.out.println("Graph potential: "+pot_graph);
//        System.out.println("Scaled: "+pot_graph/currentScale);
//        System.out.println("Tree potential:  "+pot_tree);
//        System.out.println("Scaled:  "+pot_tree/currentScale);
        return pot_graph + pot_tree;
    }
    public HashedGridFlow grad_potential(HashedGridFlow currentFlow, GridDemand demand){
//        System.out.println("Calculating Gradient in iteration "+this.iterations);
//        System.out.println("Given Flow:");
//        printFlow(currentFlow);
//        System.out.println("Given Demand:");
//        printDemand(demand);
        HashedGridFlow grad_pot_graph = currentFlow.gradient_lmax_shifted();
//        System.out.println("Gradient (Graph): ");
//        printFlow(grad_pot_graph);
//        if(printToString) System.out.println(grad_pot_graph);
//        if(printTikz) System.out.println(grad_pot_graph.tikz2D());
        // shifting nominator and denominator with same shift cancels out
        double s = -2*t.getAlpha()/t.lmax_exp_shifted_2alpha_congestion(t.getDefaultShift());
//        System.out.println("Scalar: "+s);
        t.set_edge_gradient_shift_2alpha_potential(t.getDefaultShift());
//        System.out.println("Current Tree: ");
//        printTree(t);
        GridDemand rt_times_grad = t.mult_Rt_edge_gradient();
//        System.out.println("R^T * Gradient: ");
//        printDemand(rt_times_grad);
        HashedGridFlow bt_rt_grad = rt_times_grad.toPotentialDiffEdgesFlow();
//        System.out.println("B^T * R^T * Grad: ");
//        printFlow(bt_rt_grad);
        HashedGridFlow grad_pot_tree = bt_rt_grad.scale(s);
//        System.out.println("Gradient (Tree) (-2 alpha / sum(exp(2 alpha R (b - Bf)))) * B^T * R^T * Grad:");
//        printFlow(grad_pot_tree);
        HashedGridFlow grad_pot = HashedGridFlow.add(grad_pot_graph, grad_pot_tree);
//        System.out.println("Total Gradient:");
//        printFlow(grad_pot);
        return grad_pot;
    }
    
    void printFlow(HashedGridFlow f){
        if(printToString) System.out.println(f);
        if(printTikz) {
            printTikzPre();
            System.out.println(f.tikz2D());
            printTikzPost();
        }
    }
    void printDemand(GridDemand d){
        if(printToString) System.out.println(d);
        if(printTikz) {
            printTikzPre();
            System.out.println(d.tikz2D());
            printTikzPost();
        }
    }
    void printTree(GridApproximatorTree t){
        if(printToString) System.out.println(t);
        if(printTikz) {
            printTikzPre();
            System.out.println(t.tikz2D());
            printTikzPost();
        }
    }
    
    
    
    public Tuple<HashedGridFlow, HashedGridFlow> debugAlmostRoute(GridDemand b, double eps){
        System.out.println(">>>>>>>>>> START <<<<<<<<<<<<");
        printPreamble();
        System.out.println("\\begin{document}");
        // initialization
        t.updateExcessFlows(b);
        double linf = t.linf_congestion();
        double s = ((16/eps)*Math.log(g.getN()))/(2*t.getAlpha()*linf);
        b = b.scale(s);
        currentScale = s;
        HashedGridFlow currentFlow = new HashedGridFlow(g);
        Tuple<HashedGridFlow, HashedGridFlow> iter_result;
        double delta = 0.;
        // iteration
        do{
            // [bullet 1] and pre-calculation phi, gradphi
            //System.out.println("Current Tree (before iter):");
//            printTree(t);
            System.out.println(String.format("\\section{Iteration %d with $\\delta=%f$ and $\\epsilon=%f$}", iterations, delta, eps));
            iter_result = debugiteration(currentFlow, b, eps);
            // [bullet 2] delta calculation
            delta = iter_result.b.l1();
            //System.out.println("Current scale: "+this.currentScale);
            // [bullet 3] update flow approximation
            if(delta >= eps/4){
                for(Entry<Integer, Double> e : iter_result.b.entries.entrySet()){
                    // capacity is 1 for all edges here
                    double val = (-1)*(delta/(1+4*t.getAlpha()*t.getAlpha()))*Math.signum(e.getValue());
//                    double val = -(delta/(1+4*t.getAlpha()*t.getAlpha()))*e.getValue();
                    currentFlow.add(e.getKey(), val);
                }
            }
            // just for information
            iterations++;
            // loop = "repeat" (halting condition = [bullet 4])
        } while (delta >= eps/4 && iterations < iterationLimit);
        // return
        iter_result.a = iter_result.a.scale(1./currentScale);
        iter_result.b = iter_result.b.scale(1./currentScale);
//        System.out.println("Total scaling: "+currentScale);
//        System.out.println("Stopped at iteration "+iterations);
        System.out.println(String.format("\\section{Result}"));
        System.out.println("\\subsection{Flow}");
        printFlow(iter_result.a);
        System.out.println("\\subsection{Gradient}");
        printFlow(iter_result.b);
        System.out.println("\\end{document}");
        System.out.println(">>>>>>>>>>>> END <<<<<<<<<<<<");
        return iter_result;
    }
    private Tuple<HashedGridFlow, HashedGridFlow> debugiteration(HashedGridFlow currentFlow, GridDemand demand, double eps){
        System.out.println(String.format("\\subsection{Calculation of Potential with $\\epsilon=%f$}",eps));
        System.out.println(String.format("\\subsubsection{Given Flow}"));
        printFlow(currentFlow);
        System.out.println(String.format("\\subsubsection{Given Demand}"));
        printDemand(demand);
        double pot = debugpotential(currentFlow, demand, false);
        System.out.println(String.format("Potential: %f",pot));
//        System.out.println("Pot @ "+iterations+" :   "+pot);
        while(pot < 16*(1./eps)*Math.log(g.getN())){
            currentFlow = currentFlow.scale(17./16);
            demand = demand.scale(17./16);
            this.currentScale *= 17./16;
            pot = debugpotential(currentFlow, demand, false);
        }
        System.out.println(String.format("\\subsubsection{Scaled Flow with $s=%f$}", currentScale));
        printFlow(currentFlow);
        System.out.println(String.format("\\subsubsection{Scaled Demand with $s=%f$}", currentScale));
        System.out.println(String.format("Potential after scaling with $s=%f$: %f", currentScale, pot));
        // GridFlow grad_pot = new GridFlow(g);
        HashedGridFlow grad_potential = debuggrad_potential(currentFlow, demand);
//        System.out.println("Current Flow: ");
//        printFlow(currentFlow);
//        System.out.println("Current Gradient: ");
//        printFlow(grad_potential);
//        System.out.println("");
        
        return new Tuple<>(currentFlow, grad_potential);
    }
    public double debugpotential(HashedGridFlow currentFlow, GridDemand demand, boolean print){
        if(print){
            System.out.println(String.format("Calculating Potential in iteration %d with scale %f....",this.iterations, this.currentScale));
            System.out.println("Given Flow:");
            printFlow(currentFlow);
//            System.out.println(currentFlow);
            System.out.println("Given Demand: ");
            printDemand(demand);
        }
        GridDemand bf = currentFlow.calculateExcessFlows();
        if(print){
            System.out.println("Calculated Excess Flows:");
            printDemand(bf);
        }
        GridDemand residualDemand = GridDemand.subtract(demand, bf);
        if(print){
            System.out.println("Calculated Residual Demand: ");
            printDemand(residualDemand);
        }
        double pot_graph = currentFlow.lmax_shifted();
        t.updateExcessFlows(residualDemand);
        double pot_tree = t.lmax_shifted_2alpha_congestion();
        if(print){
            System.out.println(String.format("Graph Potential: %f\\\\\n", pot_graph));
            System.out.println(String.format("Graph Potential (scaled): %f\\\\\n", pot_graph/currentScale));
            System.out.println(String.format("Tree Potential: %f\\\\\n", pot_tree));
            System.out.println(String.format("Tree Potential (scaled): %f\\\\\n", pot_tree/currentScale));
        }
        return pot_graph + pot_tree;
    }
    public HashedGridFlow debuggrad_potential(HashedGridFlow currentFlow, GridDemand demand){
        System.out.println(String.format("\\subsection{Calculating Gradient in iteration %d}",this.iterations));
        System.out.println("\\subsubsection{Given Flow}");
        printFlow(currentFlow);
        System.out.println("\\subsubsection{Given Demand}");
        printDemand(demand);
        HashedGridFlow grad_pot_graph = currentFlow.gradient_lmax_shifted();
        System.out.println("\\subsubsection{Gradient (Graph)}");
        printFlow(grad_pot_graph);
//        if(printToString) System.out.println(grad_pot_graph);
//        if(printTikz) System.out.println(grad_pot_graph.tikz2D());
        // shifting nominator and denominator with same shift cancels out
        double s = -2*t.getAlpha()/t.lmax_exp_shifted_2alpha_congestion(t.getDefaultShift());
        System.out.println("Scalar: "+s);
        t.set_edge_gradient_shift_2alpha_potential(t.getDefaultShift());
        System.out.println("\\subsubsection{Current Tree}");
        printTree(t);
        GridDemand rt_times_grad = t.mult_Rt_edge_gradient();
        System.out.println("\\subsubsection{$R^T\\cdot \\nabla \\text{lmax}$}");
        printDemand(rt_times_grad);
        HashedGridFlow bt_rt_grad = rt_times_grad.toPotentialDiffEdgesFlow();
        System.out.println("\\subsubsection{$B^T \\cdot R^T\\cdot \\nabla \\text{lmax}$}");
        printFlow(bt_rt_grad);
        HashedGridFlow grad_pot_tree = bt_rt_grad.scale(s);
        System.out.println("\\subsubsection{Gradient (Tree) $-2\\alpha\\cdot \\frac{1}{\\sum_i e^{(2\\alpha\\cdot R(b-Bf))_i}+e^{(-2\\alpha\\cdot R(b-Bf))_i}} \\cdot B^T \\cdot R^T\\cdot \\nabla \\text{lmax}$}");
        printFlow(grad_pot_tree);
        HashedGridFlow grad_pot = HashedGridFlow.add(grad_pot_graph, grad_pot_tree);
        System.out.println("\\subsubsection{Total Gradient}");
        printFlow(grad_pot);
        return grad_pot;
    }
    void printPreamble(){
        System.out.println("\\documentclass[10pt,a4paper]{article}");
        System.out.println("\\usepackage[utf8]{inputenc}");
        System.out.println("\\usepackage[english]{babel}");
        System.out.println("\\usepackage{amsmath}");
        System.out.println("\\usepackage{amsfonts}");
        System.out.println("\\usepackage{amssymb}");
        System.out.println("\\usepackage{graphicx}");
        System.out.println("\\usepackage{xcolor}");
        System.out.println("\\usepackage{tikz}");
        System.out.println("\\usetikzlibrary{positioning}");
        System.out.println("\\author{Jonas Schulz}");
    }
    void printTikzPre(){
        System.out.println("\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]");
        System.out.println("\\hspace{-4cm}");
    }
    void printTikzPost(){
        System.out.println("\\end{tikzpicture}\\\\\\\\");
    }
}
