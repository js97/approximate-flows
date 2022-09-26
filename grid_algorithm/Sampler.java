package grid_algorithm;

import java.util.Map;
import project_utils.Tuple;

/**
 *
 * @author kroka
 */
public class Sampler {
    
    static GridDemand demand1(){
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
        return d;
    }
    static GridDemand demand2(){
        GridGraph g = new GridGraph(4,4);
        GridDemand d = new GridDemand(g);
        d.set(0, +1.);
        d.set(15, -1.);
        return d;
    }
    static GridDemand demand3(){
        GridGraph g = new GridGraph(4,4);
        GridDemand d = new GridDemand(g);
        d.set(9, .5);
        d.set(6, -.4);
        d.set(5, .2);
        d.set(10, -.3);
        return d;
    }
    static GridDemand demand4(){
        GridGraph g = new GridGraph(4,4);
        GridDemand d = new GridDemand(g);
        d.set(0, +1.);
        d.set(63, -1.);
        return d;
    }
    static GridFlow flow1(){
        GridFlow f = new GridFlow(new GridGraph(4,4));
        f.set(0, 1, 1.);
        f.set(4, 8, -1.);
        f.set(3, 7, 1.);
        f.set(13, 14, -1.);
        f.set(0, 4, 1.);
        f.set(1, 5, -1.);
        return f;
    }
    public static void main(String[] args) {
        GridDemand d1 = demand3();
        GridApproximation ga = new GridApproximation(d1.g);
        GridFlow f1 = flow1();
//        HashedGridFlow s = new HashedGridFlow(f1, false);
        double eps = 0.01;
        Tuple<GridFlow, GridFlow> iter_result = ga.iteration(f1, d1, eps);
        double delta = iter_result.b.l1();
        double val = (-1)*(delta/(1+4*ga.t.getAlpha()*ga.t.getAlpha()));
        
        
        for(double h = -100; h <= 10800; h+=1){
//        for(double h = 2500; h <= 4000; h+=0.1){
            // [bullet 2] delta calculation
//            System.out.println("Current scale: "+this.currentScale);
            // [bullet 3] update flow approximation
            GridFlow f_minus_h_s = ga.flow_f_plus_h_gradsignum(iter_result.a, h*val, iter_result.b);
            
            double gh = ga.potential(f_minus_h_s, d1);
            System.out.println(h+", "+gh+", "+ga.cpg+", "+ga.cpt);
        }
    }
}
