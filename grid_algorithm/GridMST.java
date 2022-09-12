/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grid_algorithm;

import java.util.Arrays;

/**
 *
 * @author kroka
 */
public class GridMST {
    GridGraph g;
    public GridMST(GridGraph g){
        this.g = g;
    }
    public int[] iter_next(int[] current){
        return iter_next(g, current);
    }
    public GridFlow route(GridDemand b){
        return route(g, b);
    }
    public GridFlow route2(GridDemand b){
        return route2(g, b);
    }
    public static GridFlow route(GridGraph G, GridDemand b){
        GridFlow f = new GridFlow(G);
        GridDemand r = b.clone();
        int[] n = G.nodesPerDim;
        int d = n.length;
        
        for (int i = 0; i < d; i++) {
            if(n[i] == 1) continue;
            int[] v_start = new int[d];
            int[] v_end = new int[d];
            for (int j = 0; j < d; j++) {
                if(j < i){
                    v_start[j] = 0;
                    v_start[j] = 0;
                } else if(j == i){
                    v_start[j] = n[j] - 1;
                    v_end[j] = 1;
                } else {
                    v_start[j] = n[j] - 1;
                    v_end[j] = 0;
                }
            }
            boolean last_reached = false;
            for(int[] v = v_start; !last_reached; v = iter_next(G, v)){
                
                int[] u = Arrays.copyOf(v, v.length);
                u[i] = v[i]-1;
//                double V = r.get(v);
                f.set(u, v, r.get(v));
                r.set(u, r.get(v)+r.get(u));
                r.set(v, 0.);
//                System.out.println("Changed "+Arrays.toString(u)+" to "+Arrays.toString(v)+" to "+V);
                if(Arrays.equals(v, v_end)){
                    last_reached = true;
                }
            }
            
        }
        
        return f;
    }
    static int[] iter_next(GridGraph G, int[] current){
        int last_nonzero = current.length-1;
        while(last_nonzero >= 0 && current[last_nonzero] == 0){
            last_nonzero--;
        }
        int[] next = Arrays.copyOf(current, current.length);
        if(last_nonzero >= 0){
            next[last_nonzero]--;
            for (int i = last_nonzero+1; i < current.length; i++) {
                next[i] = G.nodesPerDim[i]-1;
            }
            return next;
        } else {
            return null;
        }
    }
    public static GridFlow route2(GridGraph G, GridDemand b){
        GridFlow f = new GridFlow(G);
        GridDemand r = b.clone();
        int n = G.getN();
        int i = 0;
        for (int k = n - 1; k >= 1; k--) {
            int[] v = G.toPosition(k);
            while(v[i] == 0){
                i++;
            }
            v[i]--;
//            double V = r.get(k);
            int u = G.toIndex(v);
//            System.out.println("Changed "+u+" to "+k+" to "+V);
            f.set(u, k, r.get(k));
            r.set(u, r.get(u)+r.get(k));
            r.set(k, 0.);
        }
        return f;
    }
}
