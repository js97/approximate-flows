package grid_algorithm;

import com.opencsv.CSVWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
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
    
    boolean activateLoopDetector = false;
            
    double[] loopDetector = new double[16];
    int loopIndex = 0;
    boolean loopCheck(double value){
        if(!activateLoopDetector) return false;
        for(int i = 0; i < 16; i++){
            if(16 >= iterations) break;
            if(value == loopDetector[i]){
                int loopsize = loopIndex - i;
                if(loopsize <= 0) loopsize += 16;
                System.out.println(String.format("Loop of size %d detected at iteration %d. (Value: %f)",loopsize,iterations, value));
                loopDetected = true;
                return true;
            }
        }
        loopDetector[loopIndex++] = value;
        loopIndex &= 15;
        return false;
    }
    boolean loopDetected = false;
    
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
                this.fw = new FileWriter(input_identifier);
            } catch (IOException ex) {
                Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        this.g = g;
        this.t = new GridApproximatorTree(g);
    }
    
    boolean printTikz = true;
    boolean printToString = false;
    int iterationLimit = 1_000_000;
    // debug purpose
    CSVWriter writer, iterationWriter;
    FileWriter fw;
    boolean writeToFile = false;
    
    public Tuple<GridFlow, GridFlow> AlmostRoute_inputs(GridDemand b, double eps, double alpha, DoubleSequence h){
        t.alpha = alpha;
        stepSize = h;
        return AlmostRoute(b, eps);
    }
    
    // dynamic search for optimal scaling
    boolean dynamic_opt_stepsize = true;
    double dynamic_opt_stepsize_precision = 0.001;
    // just to record how many iterations have been made
    int iterations = 0;
    public Tuple<GridFlow, GridFlow> AlmostRoute(GridDemand b, double eps){
        // initialization
        t.updateExcessFlows(b);
        double linf = t.linf_congestion();
        double s = ((16/eps)*Math.log(g.getN()))/(2*t.getAlpha()*linf);
        b = b.scale(s);
        currentScale = s;
//        try {
//            fw.write(String.format("16/eps: %f\n",(16/eps)));
//            fw.write(String.format("log(n): %f\n", Math.log(g.getN())));
//            fw.write(String.format("linf: %f\n",linf));
//            fw.write(String.format("2alpha * linf: %f\n", (2*t.getAlpha()*linf)));
//            fw.write(String.format("Term 1: %f",((16/eps)*Math.log(g.getN()))));
//            fw.write(String.format("Term 2: %f\n", (2*t.getAlpha()*linf)));
//            fw.write(String.format("Total: %f\n", ((16/eps)*Math.log(g.getN()))/(2*t.getAlpha()*linf)));
//            fw.write(String.format("Starting scale: %f\n",s));
//        } catch (IOException ex) {
//            Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
//        }
        GridFlow currentFlow = new GridFlow(g);
        Tuple<GridFlow, GridFlow> iter_result;
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
                if(dynamic_opt_stepsize)
                    val *= optimize_stepsize(iter_result.a, iter_result.b, val, b, dynamic_opt_stepsize_precision);
                else
                    val *= stepSize.at(iterations);
                for(Entry<Integer, Double> e : iter_result.b.entries.entrySet()){
                    // capacity is 1 for all edges here
//                    double val = -(delta/(1+4*t.getAlpha()*t.getAlpha()))*e.getValue();
                    currentFlow.add(e.getKey(), val * Math.signum(e.getValue()));
                }
            }
//            System.out.println("Linf residual: "+GridDemand.subtract(b, iter_result.a.calculateExcessFlows()).linf());
            if(loopCheck(iter_result.b.l1())) break;
            // just for information
            iterations++;
//            System.out.println(String.format("Delta = %f", delta));
            debug_csv : {
//                if(iterations == 0){
//                    System.out.println("Flow "+currentFlow.tikz2D());
//                    System.out.println("");
//                }
//                try {
////                    fw.write(String.format("Number of Iterations: %d\n", iterations));
////                    fw.write("STOP==========================================\n");
//                    fw.write(String.format("scale: %f  ", currentScale));
//                    fw.write(String.format("treepot: %f  ", cpt));
//                    fw.write(String.format("graphpot: %f  ", cpg));
//                    fw.write(String.format("potential: %f  ", cpt+cpg));
//                    fw.write(String.format("gnorm: %f  ", iter_result.b.l1()));
//                    fw.write(String.format("gnorm (over flow grad): %f  ", iter_result.a.gradient_lmax_shifted().l1()));
//                    fw.write(String.format("gnorm (delta): %f  ", delta));
//                    fw.write(String.format("cutbound: %f  ", iter_result.a.calculateExcessFlows().l1()));
//                    fw.write("\n");
//                    writer.writeNext(new String[]{""+delta});
//                }
//                // loop = "repeat" (halting condition = [bullet 4])
//                catch (IOException ex) {
//                    Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
//                }
            }
//            System.out.println("Iteration "+iterations+" finished.");
//            System.out.println("  Current Scale: "+currentScale);
//            System.out.println("  Current Delta: "+delta);
//            System.out.println("  Current Gradient:\n"+iter_result.b.tikz2D());
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
//                fw.write(String.format("Number of Iterations: %d\n", iterations));
//                fw.write("STOP==========================================\n");
//                fw.write(String.format("scale: %f  ", currentScale));
//                fw.write(String.format("treepot: %f  ", cpt));
//                fw.write(String.format("graphpot: %f  ", cpg));
//                fw.write(String.format("potential: %f  ", cpt+cpg));
//                fw.write(String.format("gnorm: %f  ", iter_result.b.l1()));
//                fw.write(String.format("gnorm (over flow grad): %f  ", this.grad_potential(currentFlow, b).l1()));
//                fw.write(String.format("cutbound: %f  ", iter_result.a.calculateExcessFlows().l1()));
//                fw.write(String.format("dembound: %f  ", args));
                writer.flush();
                iterationWriter.flush();
//                fw.flush();
            } catch (IOException ex) {
                Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        return iter_result;
    }
    boolean golden_section_search = true;
    double optimize_stepsize(GridFlow flow, GridFlow gradient, double std, GridDemand demand, double precision){
        double grow = 2.;
        // Initial search for interval with exponential grow
        double factor0 = 0.25;
        double factor1 = factor0 * grow;
        double factor2 = factor1 * grow;
        double pot0 = potential(flow_f_plus_h_gradsignum(flow, factor0*std, gradient), demand);
        double pot1 = potential(flow_f_plus_h_gradsignum(flow, factor1*std, gradient), demand);
        double pot2 = potential(flow_f_plus_h_gradsignum(flow, factor2*std, gradient), demand);
//        String pots_str = "Interval search: ("+factor0+", "+pot0+") - ("
//                + factor1+", "+pot1+")"+
//                + factor2+", "+pot2+")";
        while(pot1 > pot2){
            factor0 = factor1;
            factor1 = factor2;
            factor2 *= grow;
            pot0 = pot1;
            pot1 = pot2;
            pot2 = potential(flow_f_plus_h_gradsignum(flow, factor2*std, gradient), demand);
//            pots_str += " - ("+factor2+", "+pot2+")";
        }
        
        
        
        double approx_min = golden_section_search ? 
                optimize_stepsize_gss(flow, gradient, std, demand, precision, factor0, factor2) :
                optimize_stepsize_p4(flow, gradient, std, demand, precision, factor0, factor2);
//        System.out.println("Optimization step at iteration "+iterations+" with val="+std);
//        System.out.println(pots_str);
//        System.out.println(String.format("scale [%f, %f, %f, %f]",fi0,fi1,fi2,fi3));
//        System.out.println(String.format("pots [%f, %f, %f, %f]",poti0,poti1,poti2,poti3));
        return approx_min;
    }
    
    double optimize_stepsize_p4(GridFlow flow, GridFlow gradient, double std, GridDemand demand, double precision, double a, double b/*, double ya, double yb*/){
        // Refined minimum search by subsequent interval trisection
        double fi0 = a;
        double fi1 = a + (1/3.)*(b - a);
//        double fi1 = factor1;
        double fi2 = a + (2/3.)*(b - a);
        double fi3 = b;
//        double poti0 = pot0;
//        double poti1 = potential(flow_f_plus_h_grad(flow, fi1*std, gradient), demand);
//        double poti1 = pot1;
        double poti1 = potential(flow_f_plus_h_gradsignum(flow, fi1*std, gradient), demand);
        double poti2 = potential(flow_f_plus_h_gradsignum(flow, fi2*std, gradient), demand);
//        double poti3 = pot2;
        // order: fi0 < fi1 < fi2 < fi3
        while(Math.abs(poti2 - poti1) > Math.abs(precision * poti1)){
            if(poti2 > poti1){
                // minimum not at fi2, hence also not at [fi2, fi3]
                // set interval to [fi0, fi2]
                fi3 = fi2;
//                poti3 = poti2;
            } else {
                // minimum not at fi1, hence also not at [fi0, fi1]
                // set interval to [fi1, fi3]
                fi0 = fi1;
//                poti0 = poti1;
            }
            // subdivide: same formulae for both cases
            fi1 = fi0 + (1/3.)*(fi3 - fi0);
            fi2 = fi0 + (2/3.)*(fi3 - fi0);
            poti1 = potential(flow_f_plus_h_gradsignum(flow, fi1*std, gradient), demand);
            poti2 = potential(flow_f_plus_h_gradsignum(flow, fi2*std, gradient), demand);
        }
        return poti1 < poti2 ? fi1 : fi2;
    }
    double optimize_stepsize_gss(GridFlow flow, GridFlow gradient, double std, GridDemand demand, double precision, double a, double b/*, double ya, double yb*/){
        double h = b - a;
        double sq5 = Math.sqrt(5);
//        double gs = (1+sq5)/2;
//        double inv_gs = gs-1;
        double inv_gs = (sq5-1)/2;
        double inv_gs2 = (3 - sq5)/2;
        
        // search interval is [a,b]
        double x0 = a;
        double x3 = b;
//        double y0 = ya;
//        double y3 = yb;
        // for |I_(i+1)|=(1/gs)*|I_i|, we set:
        double x1 = x0 + inv_gs2 * h; // = x3 - inv_gs * h
        double x2 = x0 + inv_gs * h;
        double y1 = potential(flow_f_plus_h_gradsignum(flow, x1*std, gradient),demand);
        double y2 = potential(flow_f_plus_h_gradsignum(flow, x2*std, gradient),demand);
        
        while(Math.abs(y2-y1) > Math.abs(precision * y1)){
            if(y2 > y1){
                // minimum not in [x2, x3]
                // shrink to interval [x0, x2]
                x3 = x2;
//                y3 = y2;
                // can reuse point x1 as new x2 (golden section)
                x2 = x1;
                y2 = y1;
                // update interval size
                h = x3 - x0;
                // update point x1
                x1 = x0 + inv_gs2 * h;
                y1 = potential(flow_f_plus_h_gradsignum(flow, x1*std, gradient), demand);
            } else {
                // minimum not in [x0, x1]
                // shrink to interval [x1, x3]
                x0 = x1;
//                y0 = y1;
                // can reuse point x2 as new x1 (golden section)
                x1 = x2;
                y1 = y2;
                // update interval size
                h = x3 - x0;
                // update point x2
                x2 = x1 + inv_gs * h;
                y2 = potential(flow_f_plus_h_gradsignum(flow, x2*std, gradient), demand);
            }
        }
        
        return y1 < y2 ? x1 : x2;
    }
    GridFlow flow_f_plus_h_gradsignum(GridFlow flow, double h, GridFlow gradient){
        GridFlow cflow = new GridFlow(flow, true);
        for(Entry<Integer, Double> e : gradient.entries.entrySet()){
            // capacity is 1 for all edges here
//                    double val = -(delta/(1+4*t.getAlpha()*t.getAlpha()))*e.getValue();
            cflow.add(e.getKey(), h * Math.signum(e.getValue()));
        }
        return cflow;
    }
    
    public GridFlow CompleteRoute_inputs(GridDemand b, double eps, double alpha, DoubleSequence h){
        t.alpha = alpha;
        stepSize = h;
        return CompleteRoute(b, eps);
    }
    public GridFlow CompleteRoute(GridDemand b, double eps){
        // first route: AlmostRoute
        Tuple<GridFlow, GridFlow> f0_ = AlmostRoute(b, eps);
        int m = g.getM();
        //eventually add 1, as neither this thesis nor [She13] state whether to floor or ceil the value
        int T = (int)(Math.log(2*m)/Math.log(2));
        GridDemand bi = b;
        GridFlow fi = f0_.a;
        GridFlow f_sum = fi;
        for(int i = 1; i <= T; i++){
            // consecutive routings: AlmostRoute with residual and eps=0.5
            bi = GridDemand.subtract(bi, fi.calculateExcessFlows());
            Tuple<GridFlow, GridFlow> fi_ = AlmostRoute(bi, 0.5);
            fi = fi_.a;
            f_sum = GridFlow.add(f_sum, fi);
        }
        //TODO: Testing
        GridFlow treeRouteRes = GridMST.route2(g, bi);
        return GridFlow.add(f_sum, treeRouteRes);
    }
    
    
    boolean alternateScaling = false;
    // part of the paper's "repeat" part, namely the scaling of f and b (bullet 1),
    // and the calculation of potential and potential gradient of f (to be used for other bullet points)
    Tuple<GridFlow, GridFlow> iteration(GridFlow currentFlow, GridDemand demand, double eps){
        double pot = potential(currentFlow, demand);
//        System.out.println("Pot @ "+iterations+" :   "+pot);
        double locScale = 1.;
//        if(iterations % 50 == 0){
//            try {
//                fw.write(String.format("=======================Iteration %d =="
//                        + "==============================\n", iterations));
//            } catch (IOException ex) {
//                Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
//            }
//        }
//        System.out.println(String.format("Comparing %f with %f...",pot,16*(1./eps)*Math.log(g.getN())));
        while(pot < (alternateScaling ? 1.99 : 16) // this has to be <= 2.5125 +- 0.0125 (?)
                *(1./eps)*Math.log(g.getN())){
            // inplace scaling since caller method relies on referred objects
            double scale = alternateScaling ? 1.01 : 17./16;
            currentFlow = currentFlow.scale_inplace(scale);
            demand.scale_inplace(scale);
            this.currentScale *= scale;
            locScale *= scale;
//            try {
//                fw.append(String.format("-------------------------  Scale up! scale: %f pot: %f\n",this.currentScale, pot));
//            } catch (IOException ex) {
//                Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
//            }
            pot = potential(currentFlow, demand);
//            System.out.println("Scaling in iteration "+iterations);
//            System.out.println("  Scale is now "+currentScale);
//            System.out.println("  Potential is now "+pot);
        }
//        if(iterations == 0){
//            System.out.println(String.format("Initial Scale up was %f with potential %f", currentScale, pot));
//            System.out.println("Flow @ 0: \n"+currentFlow.tikz2D());
//            System.out.println("Demand @ 0: \n"+demand.tikz2D());
//            System.out.println("Gradient @ 0: \n"+grad_potential(currentFlow, demand).tikz2D());
//            t.set_edge_gradient_2alpha_potential();
//            System.out.println("Tree @ 0:\n"+t.tikz2D());
//        }
        if(alternateScaling){
            //Then, also scale down again (reference implementation during debug)
            while(pot > 2.1 / eps * Math.log(g.getN())){
                double scale = 0.99;
                currentFlow = currentFlow.scale_inplace(scale);
                demand.scale_inplace(scale);
                this.currentScale *= scale;
                locScale *= scale;
//                try {
//                    fw.append(String.format("-------------------------  Scale down! scale: %f pot: %f\n",this.currentScale, pot));
//                } catch (IOException ex) {
//                    Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
//                }
                pot = potential(currentFlow, demand);
            }
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
        GridFlow grad_potential = grad_potential(currentFlow, demand);
//        System.out.println("Current Flow: ");
//        printFlow(currentFlow);
//        System.out.println("Current Gradient: ");
//        printFlow(grad_potential);
        
        
        return new Tuple<>(currentFlow, grad_potential);
    }
    double cpt = 0., cpg = 0.;
    public double potential(GridFlow currentFlow, GridDemand demand){
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
    public GridFlow grad_potential(GridFlow currentFlow, GridDemand demand){
//        System.out.println("Calculating Gradient in iteration "+this.iterations);
//        System.out.println("Given Flow:");
//        printFlow(currentFlow);
//        System.out.println("Given Demand:");
//        printDemand(demand);
        /**
         * Use a respective shift for the gradient graph potential to avoid numeric error
         */
        double shift_graph = -currentFlow.linf();
        GridFlow grad_pot_graph = currentFlow.gradient_lmax_shifted(shift_graph);
        double shift = t.getDefaultShift();
//        System.out.println("Gradient (Graph): ");
//        printFlow(grad_pot_graph);
//        if(printToString) System.out.println(grad_pot_graph);
//        if(printTikz) System.out.println(grad_pot_graph.tikz2D());
        // shifting nominator and denominator with same shift cancels out
        double s = -2*t.getAlpha()/t.lmax_exp_shifted_2alpha_congestion(shift);
//        System.out.println("Scalar: "+s);
        t.set_edge_gradient_shift_2alpha_potential(shift);
//        System.out.println("Current Tree: ");
//        printTree(t);
        GridDemand rt_times_grad = t.mult_Rt_edge_gradient();
//        System.out.println("R^T * Gradient: ");
//        printDemand(rt_times_grad);
        GridFlow bt_rt_grad = rt_times_grad.toPotentialDiffEdgesFlow();
//        System.out.println("B^T * R^T * Grad: ");
//        printFlow(bt_rt_grad);
        GridFlow grad_pot_tree = bt_rt_grad.scale(s);
//        System.out.println("Gradient (Tree) (-2 alpha / sum(exp(2 alpha R (b - Bf)))) * B^T * R^T * Grad:");
//        printFlow(grad_pot_tree);
        GridFlow grad_pot = GridFlow.add(grad_pot_graph, grad_pot_tree);
//        System.out.println("Total Gradient:");
//        printFlow(grad_pot);
//        if(iterations == 1){
//            System.out.println("Gradient Calculation >");
//            System.out.println("  Input:");
//            System.out.println("    Flow:\n"+currentFlow.tikz2D());
//            System.out.println("    demand:\n"+demand.tikz2D());
//            System.out.println("    tree:\n"+t.tikz2D());
//            System.out.println("  Gradient Graph:\n"+grad_pot_graph.tikz2D());
//            System.out.println("  Scaling with "+s);
//            System.out.println("Excluded shift factor: "+t.lmax_exp_shifted_2alpha_congestion(shift));
//            System.out.println("  R*Grad(f):\n"+rt_times_grad.tikz2D());
//            System.out.println("  B*R*Grad(f):\n"+bt_rt_grad.tikz2D());
//            System.out.println("  Gradient Tree (s *B*R*Grad(f)):\n"+grad_pot_tree.tikz2D());
//            System.out.println("  Final Gradient:\n"+grad_pot.tikz2D());
//        }
        return grad_pot;
    }
    
    void printFlow(GridFlow f){
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
    
    
    
    public Tuple<GridFlow, GridFlow> debugAlmostRoute(GridDemand b, double eps){
        System.out.println(">>>>>>>>>> START <<<<<<<<<<<<");
        printPreamble();
        System.out.println("\\begin{document}");
        // initialization
        t.updateExcessFlows(b);
        double linf = t.linf_congestion();
        double s = ((16/eps)*Math.log(g.getN()))/(2*t.getAlpha()*linf);
        b = b.scale(s);
        currentScale = s;
        GridFlow currentFlow = new GridFlow(g);
        Tuple<GridFlow, GridFlow> iter_result;
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
    private Tuple<GridFlow, GridFlow> debugiteration(GridFlow currentFlow, GridDemand demand, double eps){
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
        GridFlow grad_potential = debuggrad_potential(currentFlow, demand);
//        System.out.println("Current Flow: ");
//        printFlow(currentFlow);
//        System.out.println("Current Gradient: ");
//        printFlow(grad_potential);
//        System.out.println("");
        
        return new Tuple<>(currentFlow, grad_potential);
    }
    public double debugpotential(GridFlow currentFlow, GridDemand demand, boolean print){
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
    public GridFlow debuggrad_potential(GridFlow currentFlow, GridDemand demand){
        System.out.println(String.format("\\subsection{Calculating Gradient in iteration %d}",this.iterations));
        System.out.println("\\subsubsection{Given Flow}");
        printFlow(currentFlow);
        System.out.println("\\subsubsection{Given Demand}");
        printDemand(demand);
        GridFlow grad_pot_graph = currentFlow.gradient_lmax_shifted();
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
        GridFlow bt_rt_grad = rt_times_grad.toPotentialDiffEdgesFlow();
        System.out.println("\\subsubsection{$B^T \\cdot R^T\\cdot \\nabla \\text{lmax}$}");
        printFlow(bt_rt_grad);
        GridFlow grad_pot_tree = bt_rt_grad.scale(s);
        System.out.println("\\subsubsection{Gradient (Tree) $-2\\alpha\\cdot \\frac{1}{\\sum_i e^{(2\\alpha\\cdot R(b-Bf))_i}+e^{(-2\\alpha\\cdot R(b-Bf))_i}} \\cdot B^T \\cdot R^T\\cdot \\nabla \\text{lmax}$}");
        printFlow(grad_pot_tree);
        GridFlow grad_pot = GridFlow.add(grad_pot_graph, grad_pot_tree);
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
