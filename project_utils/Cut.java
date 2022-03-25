/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package project_utils;

import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author kroka
 */
public class Cut {
    protected ArrayList<Triple<Integer, Integer, Double>> edges;
    public Cut(ArrayList<Triple<Integer, Integer, Double>> edges){
        this.edges = edges;
    }
    public Cut(Graph g, ArrayList<Integer> partition){
        this.edges = new ArrayList<>();
        for(Integer i : partition){
            for(Tuple<Integer, Double> t : g.E.get(i)){
                if(!partition.contains(t.a)){
                    edges.add(new Triple<>(i, t.a, t.b));
                }
            }
        }
    }
    public double capacity(){
        double sum = 0.;
        for(Triple<Integer, Integer, Double> t : edges){
            sum += t.c;
        }
        return sum;
    }
    public static Triple<Integer, Integer, Cut> someSTCut(Graph g){
        if(g.V.isEmpty()){
            throw new IllegalArgumentException("Can't calculate some (s,t)-cut for graph without nodes.");
        }
        if(g.V.size() < 2){
            throw new IllegalArgumentException("At least 2 nodes are needed for some (s,t)-cut.");
        }
        ArrayList<Integer> A = new ArrayList<>();
        ArrayList<Integer> B = new ArrayList<>();
        for(Integer i : g.V){
            B.add(i);
        }
        A.add(B.remove(0));
        while(A.size() < g.V.size()){
            double[] weights = new double[B.size()];
            for(int i = 0; i < B.size(); i++){
                weights[i] = weight(B.get(i), A, g);
            }
            Tuple<Integer, Double> max = Math.max(weights);
            A.add(B.remove((int)max.a));
        }
        Integer last = A.remove(A.size()-1);
        Integer second_last = A.get(A.size()-1);
        Cut c = new Cut(g, A);
        return new Triple<>(last, second_last, c);
    }
    public static double weight(int v, ArrayList<Integer> A, Graph g){
        double cur = 0.;
        ArrayList<Tuple<Integer, Double>> adjList = g.E.get(v);
        for(Tuple<Integer, Double> t : adjList){
            cur += t.b;
        }
        return cur;
    }
    public static Tuple<Double, Cut> mincut(Graph g){
        double min = Double.MAX_VALUE;
        int n = g.V.size();
        Cut c = null;
        while (n >= 2){
            Triple<Integer, Integer, Cut> somecut = someSTCut(g);
            double w = somecut.c.capacity();
            if (w < min){
                c = somecut.c;
                min = w;
            }
            //contract s and t
            g = g.contract_ext(somecut.a, somecut.b);
            n--;
        }
        return new Tuple<>(min, c);
    }
}
