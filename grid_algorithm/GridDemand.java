/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grid_algorithm;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import project_utils.Tuple;

/**
 *
 * @author kroka
 */
public class GridDemand implements Iterable<Tuple<Integer, Double>> {

    HashSet<Integer> entryIndices;
    //ArrayList<Double> entries;
    double[] entries;
    GridGraph g;
    
    public GridDemand(GridGraph g){
        this.g = g;
        //entries = new ArrayList<>(n);
        entries = new double[g.getN()];
        entryIndices = new HashSet<>();
    }
    public GridDemand(HashedGridFlow f){
        this.g = f.g;
        //entries = new ArrayList<>(n);
        entries = new double[g.getN()];
        entryIndices = new HashSet<>();
    }
    @Deprecated
    public GridDemand(int[] grid_dimensions, int n){
        //this.grid_dimensions = grid_dimensions;
        //this.n = n;
        this.g = new GridGraph(grid_dimensions);
        //entries = new ArrayList<>(n);
        entries = new double[n];
        entryIndices = new HashSet<>();
    }
    
    // currently not checking if contained to spare runtime
    public void set(Integer at, Double value){
        //entries.set(((int)at), value);
        entries[(int)at] = value;
        entryIndices.add(at);
    }
    public void add(Integer at, Double value){
        if(entryIndices.contains(at)){
            //entries.set((int)at, entries.get((int)at) + value);
            entries[(int)at] += value;
        } else {
            entries[(int)at] += value;
            entryIndices.add(at);
        }
    }
    public double get(Integer at){
        //return entries.get((int)at);
        return entries[(int)at];
    }
    public double get(int[] at_index){
        //return entries.get((int)GridGraph.toIndex(at_index, grid_dimensions));
        return entries[(int)GridGraph.toIndex(at_index, g.nodesPerDim)];
    }
    
    @Override
    public Iterator<Tuple<Integer, Double>> iterator() {
        return new Iterator<Tuple<Integer, Double>>() {
            Iterator<Integer> it = entryIndices.iterator();
            @Override
            public boolean hasNext() {
                return it.hasNext();
            }

            @Override
            public Tuple<Integer, Double> next() {
                Integer i = it.next();
                //return new Tuple<>(i, entries.get(((int)i)));
                return new Tuple<>(i, entries[(int)i]);
            }
        };
    }
    
    public double lmax_exp(){
        double sum = 0;
        int non_zeros = 0;
        for(Integer e : entryIndices){
            //double val = entries.get((int)e);
            double val = entries[(int)e];
            double exp = Math.exp(val);
            // choose following line iff Math.exp(x) is cheaper then 1/x
            // sum += exp + Math.exp(-e.c);
            // else choose this line
            sum += exp + (1/exp);
            non_zeros++;
        }
        // as this is a sparse representation, but e^0 = 1, we have to add 
        // 2 times the amount of zero-entries
        sum += 2*(g.getN()-non_zeros);
        return sum;
    }
    public double lmax(){
        return Math.log(lmax_exp());
    }
    
    public double l1(){
        double sum = 0;
        for(Integer e : entryIndices){
            double val = Math.abs(entries[(int)e]);
            sum += val;
        }
        return sum;
    }
    public double linf(){
        double max = 0;
        for(Integer e : entryIndices){
            double val = Math.abs(entries[(int)e]);
            if(val > max) max = val;
        }
        return max;
    }
    
    public static GridDemand subtract(GridDemand a, GridDemand b){
        GridDemand res = new GridDemand(a.g);
        for(Tuple<Integer, Double> e : a){
            res.set((int)e.a, e.b);
        }
        for(Tuple<Integer, Double> e : b){
            res.add((int)e.a, -1*e.b);
        }
        return res;
    }

    @Override
    protected GridDemand clone(){
        GridDemand d = new GridDemand(g);
        for(Integer i : entryIndices){
            d.entryIndices.add(i);
            d.entries[i] = entries[i];
        }
        return d;
    }
    
    public GridDemand scale(double s) {
        GridDemand d = clone();
        for(Integer i : d.entryIndices){
            d.entries[i] *= s;
        }
        return d;
    }
    public GridDemand scale_inplace(double s) {
        for(Integer i : this.entryIndices){
            this.entries[i] *= s;
        }
        return this;
    }
    
    public static double scalar_prod(GridDemand a, GridDemand b){
        double sum = 0.;
        for(Integer i : a.entryIndices){
            sum += a.entries[i] * b.entries[i];
        }
        return sum;
    }
    
//    @Deprecated
//    public GridFlow toPotentialDiffEdgesFlow(){
//        GridFlow f = new GridFlow(g);
//        for(int i = 0; i < g.getN(); i++){
//            int[][] neighbours = g.getNeighbours(g.toPosition(i));
//            for(int k = 0; k < neighbours.length; k++){
//                int other = g.toIndex(neighbours[k]);
//                double dif = this.get(i) - this.get(other);
//                f.set(i, other, dif);
//            }
//        }
//        return f;
//    }
    public HashedGridFlow toPotentialDiffEdgesFlow(){
        HashedGridFlow f = new HashedGridFlow(g);
        for(int i = 0; i < g.getN(); i++){
            int[][] neighbours = g.getNeighbours(g.toPosition(i));
            for(int k = 0; k < neighbours.length; k++){
                int other = g.toIndex(neighbours[k]);
                if(i < other){
                    double dif = this.get(i) - this.get(other);
                    f.set(i, other, -dif);
//                    if(dif != 0) {
//                        System.out.println("B^T: Setting "+i+" -> "+other+" to "+dif);
//                        System.out.println("  Other: "+get(other));
//                        System.out.println("  i: "+get(i));
//                    }
                }
            }
        }
        return f;
    }

    @Override
    public String toString() {
        return Arrays.toString(this.entries);
    }
    
    public String tikz2D(){
//        String s = "";
        DecimalFormat df = new DecimalFormat("#.000");
        String s = "\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]\n";
        //String s = "\\node (anchor) {};\n";
        for(int i = 0; i < this.g.getN(); i++){
            double scale = 3;
            s += "\\node[roundnode, minimum size = 2cm] at ("+scale*g.toPosition(i)[0]+", "+scale*g.toPosition(i)[1]+") {\\textcolor{blue!28}{"+i+":} "+(get(i) == 0.0 ? "0" : "\\textcolor{red}{"+df.format(get(i))+"}")+"};\n";
        }
        s += "\\end{tikzpicture}";
        return s;
    }
    
    
}
