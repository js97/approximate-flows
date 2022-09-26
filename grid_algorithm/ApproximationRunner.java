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
    static String filename = "";
    static void demand1(){
        g = new GridGraph(4,4);
        d = new GridDemand(g);
        d.set(0, +1.);
        d.set(1, +1.);
        d.set(2, +1.);
        d.set(3, +1.);
        d.set(12, -1.);
        d.set(13, -1.);
        d.set(14, -1.);
        d.set(15, -1.);
//        return d;
        filename = "Demand1";
    }
    static void demand2(){
        g = new GridGraph(4,4);
        d = new GridDemand(g);
        d.set(0, +1.);
        d.set(15, -1.);
//        return d;
        filename = "Demand2";
    }
    static void demand3(){
        g = new GridGraph(4,4);
        d = new GridDemand(g);
        d.set(9, .5);
        d.set(6, -.4);
        d.set(5, .2);
        d.set(10, -.3);
        filename = "Demand3";
    }
    static void demand4(){
        g = new GridGraph(8,8);
        d = new GridDemand(g);
        d.set(0, +1.);
        d.set(63, -1.);
        filename = "Single-8-8";
    }
    
    static enum mode {
        run, search, stepsize
    }
    public static void main(String[] args) {
        mode m = mode.stepsize;
//        mode m = mode.search;

        demand4();
        switch(m){
            case run:
                run();
                break;
            case search:
                search();
                break;
            case stepsize:
                search_opt_step_v2();
                break;
        }
        
        
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
            GridFlow res = appr.AlmostRoute_inputs(d, 0.1, 3., hseq).a;
            System.out.println(String.format("Run took %d iterations.",appr.iterations));
            System.out.println(res.tikz2D());
            System.out.println("Demands:");
            System.out.println(res.calculateExcessFlows().tikz2D());
        }
    }
    
    static void search(){
        for(double eps = 0.01; eps <= 0.1; eps += 0.01){
            try {
                CSVWriter writer = new CSVWriter(new FileWriter("LogsFixed/"+filename+"/iterations_h_alpha_metrics2_eps_"+eps+".csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);

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
                        Tuple<GridFlow, GridFlow> result = appr.AlmostRoute_inputs(d, eps, alpha, hseq);

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
    static void search_optimize_stepsize(){
        GridApproximation appr = new GridApproximation(g, "mytreeoutput.txt");
        appr.printTikz = false;
        appr.iterationLimit = 1_000_000;
        DoubleSequence hseq = new DoubleSequence(){
                        public double at(int index){
                            return 1.;
                        }
                    };
        appr.dynamic_opt_stepsize = false;
        System.out.println("Performing standard stepsize run...");
        long timeStart = System.nanoTime();
        GridFlow res = appr.AlmostRoute_inputs(d, 0.1, 3., hseq).a;
        long timeStop = System.nanoTime();
        System.out.println(String.format("Run took %d iterations and %d time.",appr.iterations, timeStop - timeStart));
        System.out.println(res.tikz2D());
        System.out.println("Demands:");
        System.out.println(res.calculateExcessFlows().tikz2D());
        String str = "Iterations:\n  Standard: "+appr.iterations+" It. in time: "+(timeStop-timeStart);
        for(double prec = 1.001; prec >= 0.0009999; prec -= 0.005){
            GridApproximation apprOpt = new GridApproximation(g);
            apprOpt.printTikz = false;
            apprOpt.iterationLimit = 1_000_000;
            DoubleSequence hseq2 = new DoubleSequence(){
                public double at(int index){
                    return 1.;
                }
            };
            apprOpt.dynamic_opt_stepsize = true;
            apprOpt.dynamic_opt_stepsize_precision = prec;
            System.out.println("Performing dynamic stepsize optimization with precision "+prec+"...");
            timeStart = System.nanoTime();
            GridFlow res2 = apprOpt.AlmostRoute_inputs(d, 0.1, 3., hseq).a;
            timeStop = System.nanoTime();
            System.out.println(String.format("Run took %d iterations and %d time.",apprOpt.iterations, timeStop - timeStart));
            System.out.println(res2.tikz2D());
            System.out.println("Demands:");
            System.out.println(res2.calculateExcessFlows().tikz2D());
            str += "\n  P "+prec+": "+apprOpt.iterations + " It. in time: "+(timeStop-timeStart);
        }
        System.out.println(str);
    }
        static void search_opt_step_v2(){
        for(double eps = 0.01; eps <= 0.10001; eps += 0.01){
            try {
                CSVWriter writer = new CSVWriter(new FileWriter("Stepsize-Optimization/"+filename+"/iterations_time_eps_"+eps+"-p4-gss-reference.csv"), CSVWriter.DEFAULT_SEPARATOR, CSVWriter.NO_QUOTE_CHARACTER);

//                writer.writeNext(new String[]{filename});
                writer.writeNext(new String[]{"Precision", "Alpha",
                    "P4 iter", "GSS iter", "ref iter",
                    "P4 time", "GSS time", "ref time",
                    "P4 residual cong", "GSS residual cong", "ref residual cong"
                });

                int iterationLimit = 500_000;
                
                double startalpha = 0.5;
                double endalpha = 20.5;
                int alphasteps = 40;
                double small_numeric = 0.000000001;
                for(double alpha = startalpha; alpha <= endalpha + small_numeric; alpha += (endalpha-startalpha)/(alphasteps-1)){
                    GridApproximation reference = new GridApproximation(g);
                    reference.iterationLimit = iterationLimit;
                    reference.printTikz = false;
                    DoubleSequence seq1 = (i) -> 1.;
                    reference.dynamic_opt_stepsize = false;
                    System.out.println(String.format("Running reference with eps=%f, alpha=%f",eps, alpha));
                    long refTimeStart = System.nanoTime();
                    Tuple<GridFlow, GridFlow> ref_result = reference.AlmostRoute_inputs(d, eps, alpha, seq1);
                    long refTimeStop = System.nanoTime();
                    long duration_ref = refTimeStop - refTimeStart;
                    GridDemand refResidual = GridDemand.subtract(d, ref_result.a.calculateExcessFlows());
                    double refmaxentry = refResidual.linf();
//                    double reflc = reference.t.compute_lemma_congestion(d);
//                    double reflinf_graph = ref_result.a.linf();
//                    GridApproximatorTree reftreeFinal = new GridApproximatorTree(g);
//                    reftreeFinal.updateExcessFlows(refResidual);
//                    double reflinf_tree = 2*alpha*reftreeFinal.linf_congestion();
//                    writer.writeNext(new String[]{"OFF", ""+alpha, ""+reference.iterations, ""+ref_result.a.get(4, 8), ""+reflc, 
//                        ""+reference.t.compute_lemma_congestion_upper(reflc, eps), ""+refmaxentry,
//                        ""+reflinf_graph, ""+reflinf_tree, ""+(refTimeStop - refTimeStart)});

                    for(double h = 1.024; h >= 0.000012499; h/= 2){
                        GridApproximation appr = new GridApproximation(g);
                        GridApproximation appr_gss = new GridApproximation(g);
                        appr.iterationLimit = iterationLimit;
                        appr_gss.iterationLimit = iterationLimit;
                        appr.printTikz = false;
                        appr_gss.printTikz = false;
                        DoubleSequence hseq = new DoubleSequence(){
                            public double at(int index){
                                return 1;
                            }
                        };
                        appr.dynamic_opt_stepsize = true;
                        appr_gss.dynamic_opt_stepsize = true;
                        appr.dynamic_opt_stepsize_precision = h;
                        appr_gss.dynamic_opt_stepsize_precision = h;
                        appr.golden_section_search = false;
                        appr_gss.golden_section_search = true;
                        System.out.println(String.format("Running P4 with eps=%f, alpha=%f, precision=%f",eps,alpha,h));
                        long timeStart = System.nanoTime();
                        Tuple<GridFlow, GridFlow> result = appr.AlmostRoute_inputs(d, eps, alpha, hseq);
                        long timeStop = System.nanoTime();
                        long duration_p4 = timeStop - timeStart;
                        
                        GridDemand finalResidual = GridDemand.subtract(d, result.a.calculateExcessFlows());
                        double maxentry = finalResidual.linf();
//                        double lc = appr.t.compute_lemma_congestion(d);
//                        double linf_graph = result.a.linf();
//                        GridApproximatorTree treeFinal = new GridApproximatorTree(g);
//                        treeFinal.updateExcessFlows(finalResidual);
//                        double linf_tree = 2*alpha*treeFinal.linf_congestion();
//                        writer.writeNext(new String[]{""+h, ""+alpha, ""+appr.iterations, ""+result.a.get(4, 8), ""+lc, 
//                            ""+appr.t.compute_lemma_congestion_upper(lc, eps), ""+maxentry,
//                            ""+linf_graph, ""+linf_tree, ""+(timeStop - timeStart),
//                            // reference data:
//                            ""+reference.iterations, ""+(refTimeStop - refTimeStart), ""+reflc, ""+refmaxentry, ""+reflinf_graph, ""+reflinf_tree
//                        });

                        System.out.println(String.format("Running GSS with eps=%f, alpha=%f, precision=%f",eps,alpha,h));
                        timeStart = System.nanoTime();
                        Tuple<GridFlow, GridFlow> result_gss = appr_gss.AlmostRoute_inputs(d, eps, alpha, hseq);
                        timeStop = System.nanoTime();
                        double duration_gss = timeStop - timeStart;
                        GridDemand finalResidual_gss = GridDemand.subtract(d, result_gss.a.calculateExcessFlows());
                        double maxentry_gss = finalResidual_gss.linf();
                        
                        writer.writeNext(new String[]{""+h,""+alpha,
                            ""+appr.iterations, ""+appr_gss.iterations, ""+reference.iterations,
                            ""+duration_p4, ""+duration_gss, ""+duration_ref,
                            ""+maxentry, ""+maxentry_gss, ""+refmaxentry
                        });
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
}
