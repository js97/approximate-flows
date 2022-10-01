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
 * This class implements the algorithmic functionalities of the algorithm from
 * [She13], our further adjustments, and some debug functionalities.
 * @author Jonas Schulz
 */
public class GridApproximation {
    /**
     * Grid graph to operate on.
     */
    GridGraph g;
    /**
     * Approximator for the grid graph {@link #g}.
     */
    GridApproximatorTree t;
    /**
     * The current scale of the demand and flow.
     */
    double currentScale = 1.;
    /**
     * Step size to be used, depending on the iteration.
     * The {@link DoubleSequence} interface offers a simple possibility for dynamic
     * step sizes; its function is the implementation of some map from integers
     * (i.e. the domain of iteration counts) to doubles (i.e. the step size).
     */
    DoubleSequence stepSize = DoubleSequence.one;
    /**
     * Flag to (de-)activate the primitive {@link #loopDetector}.
     */
    boolean activateLoopDetector = false;
    /**
     * Primitive loop detector for debug purpose.
     * It will store 16 values to be checked for a loop.
     * For use as a ring buffer, we additionally introduce the pointer {@link #loopIndex}.
     */
    double[] loopDetector = new double[16];
    /**
     * The current index of the ring buffer {@link #loopDetector}.
     */
    int loopIndex = 0;
    /**
     * Feedback from the loop detector. <code>true</code> iff a detected loop has been registered
     * during {@link #loopCheck(double)}.
     */
    boolean loopDetected = false;
    
    /**
     * Uses the primitive loop detector iff {@link #activateLoopDetector} is <code>true</code>.
     * For this, the ring buffer given by {@link #loopDetector} and {@link #loopIndex} is used.
     * If there were more than 16 iterations and some value in the buffer has been repeated,
     * {@link #loopDetected} is set to <code>true</code> and the size of the detected loop is printed on the standard output stream, {@link System#out}.
     * @param value The value to be added into the ring buffer.
     * @return <code>true</code> iff a loop has been detected.
     */
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
    
    /**
     * Path to the folder where debug information should be printed in.
     */
    String path = "./Logs/";
    
    /**
     * Standard constructor.
     * @param g The grid graph to operate on.
     */
    public GridApproximation(GridGraph g){
        // <editor-fold defaultstate="collapsed" desc="CSV Debug Initialization">
        debug_csv : {
            String subfolder = "";
//            String path = "C:/Users/kroka/Documents/Masterarbeit/Logs/";
            try {
                this.writer = new CSVWriter(new FileWriter(path + subfolder + "deltas.csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);
                writer.writeNext(new String[]{"delta"});
                this.iterationWriter = new CSVWriter(new FileWriter(path + subfolder + "potentials.csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);
                iterationWriter.writeNext(new String[]{"potential", "current scale", String.format("edgeflow-%d-%d",4,8), "edgeflow-rescaled"});
                } catch (IOException ex) {
                Logger.getLogger(GridApproximation.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        // </editor-fold>
        this.g = g;
        this.t = new GridApproximatorTree(g);
    }
    
    /**
     * Constructor with debug purpose.
     * @param g The grid graph to operate on.
     * @param input_identifier A {@link String} containing information unique for the given input. It will be 
     * part of the output file's name, and equal input identifiers will lead to the output files being overwritten.
     */
    public GridApproximation(GridGraph g, String input_identifier){
        // <editor-fold defaultstate="collapsed" desc="CSV Debug Initialization">
        debug_csv : {
            String subfolder = "";
//            String path = "C:/Users/kroka/Documents/Masterarbeit/Logs/";
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
        // </editor-fold>
        this.g = g;
        this.t = new GridApproximatorTree(g);
    }
    
    /**
     * Flag for activating any TikZ output.
     */
    boolean printTikz = true;
    /**
     * Flag for activating any data structure string representation output.
     */
    boolean printToString = false;
    /**
     * Limit for the iterations. The program will terminate after this number of iterations,
     * regardless of whether a satisfying solution has been found.
     */
    int iterationLimit = 1_000_000;
    // debug purpose
    /**
     * Writer for debug purpose.
     */
    CSVWriter writer, iterationWriter;
    /**
     * Writer for debug purpose.
     */
    FileWriter fw;
    /**
     * Flag for writing to {@link #fw}.
     * @deprecated Old debug feature. Currently not used.
     */
    @Deprecated
    boolean writeToFile = false;
    
    /**
     * Same as {@link #AlmostRoute(GridDemand, double)}, but with additional inputs for <code>alpha</code> and <code>h</code>.
     * @param b The demand to be routed.
     * @param eps The relative precision &epsilon;.
     * @param alpha Alpha parameter for the &alpha;-approximator {@link #t}.
     * @param h Simple dynamic map from iterations to step sizes.
     * @return 
     */
    public Tuple<GridFlow, GridFlow> AlmostRoute_inputs(GridDemand b, double eps, double alpha, DoubleSequence h){
        t.alpha = alpha;
        stepSize = h;
        return AlmostRoute(b, eps);
    }
    
    /**
     * (De-)activates the dynamic step size optimization with line search.
     */
    boolean dynamic_opt_stepsize = true;
    /**
     * Sets the relative precision of the dynamic step size optimization.
     */
    double dynamic_opt_stepsize_precision = 0.001;
    /**
     * Current iteration level.
     */
    int iterations = 0;
    
    /**
     * Routes most part of b with relative precision &epsilon;.
     * The algorithm terminates prematurely if the iteration limit is reached, or
     * the loop detector is activated and detects a loop.
     * @param b The demand to be routed.
     * @param eps The relative precision &epsilon;.
     * @return A tuple of flow <i>f</i> and the gradient <i>&nabla;&phi;(f)</i>.
     */
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
    
    /**
     * Switches between golden section search (<code>true></code>) and ternary search
     * (<code>false</code>) in the dynamic step size optimization line search.
     */
    boolean golden_section_search = true;
    
    /**
     * Performs a line search to find the step size with the steepest gradient descent.
     * This method first determines the search interval, and then calls either
     * the golden section search ({@link #optimize_stepsize_gss(GridFlow, GridFlow, double, GridDemand, double, double, double)})
     * or the ternary search ({@link #optimize_stepsize_p4(GridFlow, GridFlow, double, GridDemand, double, double, double)}),
     * depending on the {@link #golden_section_search} flag.
     * @param flow Current flow.
     * @param gradient Current gradient.
     * @param std Standard step size.
     * @param demand Demanded excess flows.
     * @param precision Relative precision of the line search.
     * @return The factor to multiply the standard step size with to obtain the optimum suggested by the line search.
     */
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
    
    /**
     * Implements the ternary search.
     * @param flow Current flow.
     * @param gradient Current gradient.
     * @param std Standard step size.
     * @param demand Demanded excess flows.
     * @param precision Relative precision of the line search.
     * @param a Lower bound of the search interval.
     * @param b Higher bound of the search interval.
     * @return The factor to multiply the standard step size with to obtain the optimum suggested by the line search.
     */
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
    
    /**
     * Implements the golden section search.
     * @param flow Current flow.
     * @param gradient Current gradient.
     * @param std Standard step size.
     * @param demand Demanded excess flows.
     * @param precision Relative precision of the line search.
     * @param a Lower bound of the search interval.
     * @param b Higher bound of the search interval.
     * @return The factor to multiply the standard step size with to obtain the optimum suggested by the line search.
     */
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
    
    /**
     * Helper function to calculate the value of <i>f<sub>e</sub> + h &#8729; sgn(&nabla;&phi;<sub>e</sub>)</i> for all edges <i>e</i>.
     * @param flow <i>f</i>, the flow.
     * @param h <i>h</i>, the step size.
     * @param gradient <i>&nabla;&phi;</i>, the potential gradient.
     * @return The flow with values <i>f<sub>e</sub> + h &#8729; sgn(&nabla;&phi;<sub>e</sub>)</i> at edges <i>e</i>.
     */
    GridFlow flow_f_plus_h_gradsignum(GridFlow flow, double h, GridFlow gradient){
        GridFlow cflow = new GridFlow(flow, true);
        for(Entry<Integer, Double> e : gradient.entries.entrySet()){
            // capacity is 1 for all edges here
//                    double val = -(delta/(1+4*t.getAlpha()*t.getAlpha()))*e.getValue();
            cflow.add(e.getKey(), h * Math.signum(e.getValue()));
        }
        return cflow;
    }
    
    /**
     * Same as {@link #CompleteRoute(GridDemand, double)}, but with additional parameters for &alpha; and the step size.
     * @param b The demand to be routed.
     * @param eps The relative precision &epsilon;.
     * @param alpha Alpha parameter for the &alpha;-approximator {@link #t}.
     * @param h Simple dynamic map from iterations to step sizes.
     * @return Flow that completely routes <code>b</code>, optimal within relative precision &epsilon;.
     */
    public GridFlow CompleteRoute_inputs(GridDemand b, double eps, double alpha, DoubleSequence h){
        t.alpha = alpha;
        stepSize = h;
        return CompleteRoute(b, eps);
    }
    /**
     * Completely route a demand through a grid graph, with the approximative algorithm described in [She13].
     * @param b The demand to be routed.
     * @param eps The relative precision &epsilon;.
     * @return Flow that completely routes <code>b</code>, optimal within relative precision &epsilon;.
     */
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
    
    /**
     * Alternate the scaling for test purposes.
     */
    boolean alternateScaling = false;
    /**
     * Main part of the iteration workload.
     * Performs the scaling and calculates the gradient.
     * @param currentFlow The current flow.
     * @param demand The demand to be routed.
     * @param eps The relative precision the final solution should have.
     * @return The scaled flow, together with the calculated gradient <i>&nabla;&phi;(f)</i>.
     */
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
    
    /**
     * Current tree potential component <i>lmax(2&alpha;R(b-Bf))</i>.
     */
    double cpt = 0.;
    /**
     * Current graph potential component <i>lmax(f)</i>.
     */
    double cpg = 0.;
    
    /**
     * Calculates <i>&phi;(f) = lmax(f) + lmax(2&alpha;R(b-Bf))</i>.
     * @param currentFlow <i>f</i>.
     * @param demand <i>b</i>.
     * @return <i>&phi;(f)</i>.
     */
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
    
    /**
     * Calculates <i>&nabla;&phi;(f) = &nabla;lmax(f) - 2&alpha;B<sup>T</sup>R<sup>T</sup>&nabla;lmax(2&alpha;R(b-Bf))</i>.
     * @param currentFlow <i>f</i>.
     * @param demand <i>b</i>.
     * @return <i>&nabla;&phi;(f)</i>.
     */
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
    
    /**
     * Prints String/TikZ representations of a flow on {@link System#out}, according to the flags {@link #printToString} and {@link #printTikz}.
     * @param f The flow to be printed.
     */
    void printFlow(GridFlow f){
        if(printToString) System.out.println(f);
        if(printTikz) {
            printTikzPre();
            System.out.println(f.tikz2D());
            printTikzPost();
        }
    }
    
    /**
     * Prints String/TikZ representations of a demand on {@link System#out}, according to the flags {@link #printToString} and {@link #printTikz}.
     * @param d The demand to be printed.
     */
    void printDemand(GridDemand d){
        if(printToString) System.out.println(d);
        if(printTikz) {
            printTikzPre();
            System.out.println(d.tikz2D());
            printTikzPost();
        }
    }
    
    /**
     * Prints String/TikZ representations of an approximator tree on {@link System#out}, according to the flags {@link #printToString} and {@link #printTikz}.
     * @param f The approximator tree to be printed.
     */
    void printTree(GridApproximatorTree t){
        if(printToString) System.out.println(t);
        if(printTikz) {
            printTikzPre();
            System.out.println(t.tikz2D());
            printTikzPost();
        }
    }
    
    /**
     * @deprecated Old debug functionality.
     */
    @Deprecated
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
    /**
     * @deprecated Old debug functionality.
     */
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
    /**
     * @deprecated Old debug functionality.
     */
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
    /**
     * @deprecated Old debug functionality.
     */
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
    
    /**
     * Prints the necessary imports for LaTeX to enable the functionality
     * of the generated TikZ-representations. You still need to begin and end
     * the document with <code>\begin{document}</code> and <code>\end{document}</code>.
     */
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
    /**
     * Prints the Prefix for the TikZ-picture.
     */
    void printTikzPre(){
        System.out.println("\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]");
        System.out.println("\\hspace{-4cm}");
    }
    /**
     * Prints the Suffix for the TikZ-picture.
     */
    void printTikzPost(){
        System.out.println("\\end{tikzpicture}\\\\\\\\");
    }
}
