/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grid_algorithm;

import java.util.Map.Entry;
import project_utils.Tuple;

/**
 *
 * @author kroka
 */
public class GridApproximation {
    GridGraph g;
    GridApproximatorTree t;
    double currentScale = 1.;
    public GridApproximation(GridGraph g){
        this.g = g;
        this.t = new GridApproximatorTree(g);
    }
    
    // just to record how many iterations have been made
    private int iterations = 0;
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
            iter_result = iteration(currentFlow, b, eps);
            // [bullet 2] delta calculation
            delta = iter_result.b.l1();
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
        } while (delta >= eps/4 && iterations < 1500);
        // return
        iter_result.a = iter_result.a.scale(1./currentScale);
        iter_result.b = iter_result.b.scale(1./currentScale);
        System.out.println("Total scaling: "+currentScale);
        System.out.println("Stopped at iteration "+iterations);
        return iter_result;
    }
    // part of the paper's "repeat" part, namely the scaling of f and b (bullet 1),
    // and the calculation of potential and potential gradient of f (to be used for other bullet points)
    private Tuple<HashedGridFlow, HashedGridFlow> iteration(HashedGridFlow currentFlow, GridDemand demand, double eps){
        double pot = potential(currentFlow, demand);
        System.out.println("Pot @ "+iterations+" :   "+pot);
        while(pot < 16*(1./eps)*Math.log(g.getN())){
            currentFlow = currentFlow.scale(17./16);
            demand = demand.scale(17./16);
            this.currentScale *= 17./16;
            pot = potential(currentFlow, demand);
        }
        // GridFlow grad_pot = new GridFlow(g);
        HashedGridFlow grad_potential = grad_potential(currentFlow, demand);
        System.out.println("Current Flow: "+currentFlow);
        System.out.println("Current Gradient: "+grad_potential);
        
        
        return new Tuple<>(currentFlow, grad_potential);
    }
    public double potential(HashedGridFlow currentFlow, GridDemand demand){
        GridDemand bf = currentFlow.calculateExcessFlows();
        GridDemand residualDemand = GridDemand.subtract(demand, bf);
        double pot_graph = currentFlow.lmax_shifted();
        t.updateExcessFlows(residualDemand);
        double pot_tree = t.lmax_shifted_2alpha_congestion();
        System.out.println("Graph potential: "+pot_graph);
        System.out.println("Tree potential:  "+pot_tree);
        return pot_graph + pot_tree;
    }
    public HashedGridFlow grad_potential(HashedGridFlow currentFlow, GridDemand demand){
        HashedGridFlow grad_pot_graph = currentFlow.gradient_lmax_shifted();
        // shifting nominator and denominator with same shift cancels out
        double s = -2*t.getAlpha()/t.lmax_exp_shifted_2alpha_congestion(t.getDefaultShift());
        t.set_edge_gradient_shift_2alpha_potential(t.getDefaultShift());
        GridDemand rt_times_grad = t.mult_Rt_edge_gradient();
        HashedGridFlow bt_rt_grad = rt_times_grad.toPotentialDiffEdgesFlow();
        HashedGridFlow grad_pot_tree = bt_rt_grad.scale(s);
        HashedGridFlow grad_pot = HashedGridFlow.add(grad_pot_graph, grad_pot_tree);
        return grad_pot;
    }
}
