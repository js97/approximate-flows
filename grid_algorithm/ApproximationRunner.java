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
    static String tikz = "\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]\n"
            + "%s\n\\end{tikzpicture}";
    public static void main(String[] args) {
        GridGraph g = new GridGraph(4,4);
        GridDemand d = new GridDemand(g);
        d.set(0, +1.);
        d.set(1, +1.);
        d.set(2, +1.);
        d.set(3, +1.);
        d.set(12, -1.);
        d.set(13, -1.);
        d.set(14, -1.);
        d.set(15, -1.);
        
        try {
            CSVWriter writer = new CSVWriter(new FileWriter("Demand2/iterations_per_h_alpha_combined_eps_01.csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);
            
            int iterationLimit = 100_000;
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
                    double eps = 0.05;
                    Tuple<HashedGridFlow, HashedGridFlow> result = appr.AlmostRoute_inputs(d, eps, alpha, hseq);
                    
                    writer.writeNext(new String[]{""+h, ""+alpha, ""+appr.iterations, ""+result.a.get(4, 8)});
                    
                    
//                    if(alpha < 8) alpha += 0.25;
//                    else alpha += 1;
                }
            }
            writer.flush();
        } catch (IOException ex) {
            Logger.getLogger(ApproximationRunner.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        System.exit(0);
    }
}
