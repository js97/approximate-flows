/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grid_algorithm;

import com.opencsv.CSVWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import project_utils.DoubleSequence;
import project_utils.Tuple;

/**
 *
 * @author kroka
 */
public class ApproximationRunner {
    static GridGraph g;
    static GridDemand d;
    static String tikz = "\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]\n"
            + "%s\n\\end{tikzpicture}";
    public static void main(String[] args) {
        g = new GridGraph(4,4);
        d = new GridDemand(g);
        d.set(0, +1.);
//        d.set(1, +1.);
//        d.set(2, +1.);
//        d.set(3, +1.);
//        d.set(12, -1.);
//        d.set(13, -1.);
//        d.set(14, -1.);
        d.set(15, -1.);
        
        
//        search();
        run();
        
        
        System.exit(0);
    }
    
    static void run(){
        {
//            GridApproximation appr = new GridApproximation(g);
            GridApproximation appr = new GridApproximation(g, "mytreeoutput.txt");
            appr.printTikz = false;
            appr.iterationLimit = 1_000_000;
            DoubleSequence hseq = new DoubleSequence(){
                            public double at(int index){
                                return 1.;
                            }
                        };
            HashedGridFlow res = appr.AlmostRoute_inputs(d, 0.01, 3., hseq).a;
            System.out.println(String.format("Run took %d iterations.",appr.iterations));
            System.out.println(res.tikz2D());
            System.out.println("Demands:");
            System.out.println(res.calculateExcessFlows().tikz2D());
        }
    }
    
    static void search(){
        try {
            double eps = 0.1;
            CSVWriter writer = new CSVWriter(new FileWriter("Demand2/iterations_h_alpha_metrics2_loopdet__eps_01.csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);
            
            int iterationLimit = 500_000;
            for(double h = 0.4; h <= 1.3; h+= 0.05){
                double startalpha = 0.15;
                double endalpha = 20;
                int alphasteps = 34;
                double small_numeric = 0.000000001;
                for(double alpha = startalpha; alpha <= endalpha + small_numeric; alpha += (endalpha-startalpha)/(alphasteps-1)){
                    GridApproximation appr = new GridApproximation(g);
                    appr.iterationLimit = iterationLimit;
                    appr.printTikz = false;
                    final double h2 = h;
                    DoubleSequence hseq = new DoubleSequence(){
                        public double at(int index){
                            return h2;
                        }
                    };
                    System.out.println(String.format("Running with eps=%f, alpha=%f, h=%f",eps,alpha,h2));
                    Tuple<HashedGridFlow, HashedGridFlow> result = appr.AlmostRoute_inputs(d, eps, alpha, hseq);
                    
                    GridDemand finalResidual = GridDemand.subtract(d, result.a.calculateExcessFlows());
                    double maxentry = finalResidual.linf();
                    double lc = appr.t.compute_lemma_congestion(d);
                    double linf_graph = result.a.linf();
                    GridApproximatorTree treeFinal = new GridApproximatorTree(g);
                    treeFinal.updateExcessFlows(finalResidual);
                    double linf_tree = 2*alpha*treeFinal.linf_congestion();
                    writer.writeNext(new String[]{""+h, ""+alpha, ""+appr.iterations, ""+result.a.get(4, 8), ""+lc, 
                        ""+appr.t.compute_lemma_congestion_upper(lc, eps), ""+maxentry, ""+appr.loopDetected,
                        ""+linf_graph, ""+linf_tree});
                    
                    
//                    if(alpha < 8) alpha += 0.25;
//                    else alpha += 1;
                }
            }
            writer.flush();
        } catch (IOException ex) {
            Logger.getLogger(ApproximationRunner.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
