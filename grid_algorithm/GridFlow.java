/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grid_algorithm;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;
import project_utils.Triple;
import project_utils.Tuple;

/**
 *
 * @author kroka
 */
public class GridFlow implements Iterable<Triple<Integer,Integer,Double>>{
    HashSet<Triple<Integer, Integer, Double>> entries = new HashSet<>();
    int[] grid_dimensions;
    int n;
    int m;
    
    public GridFlow(GridGraph g){
        this.grid_dimensions = g.nodesPerDim;
        this.n = g.getN();
        calculateM();
    }
    private GridFlow(GridFlow f, boolean clone){
        this.entries = new HashSet<>();
        this.grid_dimensions = f.grid_dimensions;
        this.n = f.n;
        this.m = f.m;
        if(clone){
            for(Triple<Integer, Integer, Double> t : f.entries){
                this.entries.add(new Triple<>((int)t.a, (int)t.b, (double)t.c));
            }
        }
    }
    
    public int calculateM(){
        int sum = 0;
        for(int i = 0; i < grid_dimensions.length; i++){
            sum += (n/grid_dimensions[i])*(grid_dimensions[i]-1);
        }
        return 2*sum;
    }

    public int getM() {
        return m;
    }

    public int getN() {
        return n;
    }
    

    @Override
    public Iterator iterator() {
        return entries.iterator();
    }
    public void set(Integer from, Integer to, Double value){
        entries.add(new Triple<>(from, to, value));
    }
    public void set(int[] from, int[] to, Double value){
        int from_int = GridGraph.toIndex(from, this.grid_dimensions);
        int to_int = GridGraph.toIndex(to, this.grid_dimensions);
        this.set(from_int, to_int, value);
    }
    
    // this method represents the calculation B*f for 
    // - B the incidence matrix
    // - f the flow
    public GridDemand calculateExcessFlows(){
        GridDemand b = new GridDemand(this);
        for(Triple<Integer, Integer, Double> e : entries){
            b.add(e.a, -e.c);
            b.add(e.b, e.c);
        }
        return b;
    }
    
    public double lmax_exp(){
        double sum = 0;
        int non_zeros = 0;
        for(Triple<Integer, Integer, Double> e : entries){
            double exp = Math.exp(e.c);
            // choose following line iff Math.exp(x) is cheaper then 1/x
            // sum += exp + Math.exp(-e.c);
            // else choose this line
            sum += exp + (1/exp);
            non_zeros++;
        }
        // as this is a sparse representation, but e^0 = 1, we have to add 
        // 2 times the amount of zero-entries
        sum += 2*(m-non_zeros);
        return sum;
    }
    public double lmax(){
        return Math.log(lmax_exp());
    }
    public GridFlow gradient_lmax(){
        double sum_exps = lmax_exp();
        GridFlow f = new GridFlow(this, false);
        for(Triple<Integer, Integer, Double> t : this.entries){
            double exp = Math.exp(t.c);
            f.set((int)t.a, (int)t.b, (exp-(1./exp))/sum_exps);
        }
        // where there is i with no entry, xi = 0, hence e^xi - e^(-xi) = 0, thus 
        // also the gradient is 0 at i
        return f;
    }

    public GridFlow scale(double d) {
        GridFlow cl = new GridFlow(this, true);
        for(Triple<Integer, Integer, Double> t : cl.entries){
            t.c = t.c * d;
        }
        return cl;
    }
    public static GridFlow add(GridFlow a, GridFlow b) {
        GridFlow cl = new GridFlow(a, true);
        for(Triple<Integer, Integer, Double> t : a.entries){
            cl.set(t.a, t.b, t.c);
        }
        for(Triple<Integer, Integer, Double> t : b.entries){
            cl.set(t.a, t.b, t.c);
        }
        cl.contractDuplicates(true);
        return cl;
    }
    void contractDuplicates(boolean contractBidirectional){
        HashMap<Integer, Double> hm = new HashMap(m);
        for(Triple<Integer, Integer, Double> t : entries){
            Triple<Integer, Integer, Double> tohash = (contractBidirectional && t.a > t.b) ? new Triple<>(t.b, t.a, -t.c) : t;
            int hash = hashEntry(tohash);
            if(!hm.containsKey(hash))
                hm.put(hash, tohash.c);
            else {
                hm.put(hash, hm.get(hash) + tohash.c);
            }
        }
        entries = new HashSet<>();
        for(Entry<Integer, Double> e : hm.entrySet()){
            if(e.getValue() != 0){
                Tuple<Integer, Integer> pos = hashToTuple(e.getKey());
                entries.add(new Triple<>(pos.a, pos.b, e.getValue()));
            }
        }
    }
    
    int hashEntry(Triple<Integer, Integer, Double> t){
        return (((int)t.a) << 16) + t.b;
    }
    Tuple<Integer, Integer> hashToTuple(int hash){
        int second = hash & 0x0000_ffff;
        int first = (hash >> 16) & 0x0000_ffff;
        return new Tuple<>(first, second);
    }

    public double l1() {
        double curr = 0;
        for(Triple<Integer, Integer, Double> e : entries){
            curr += Math.abs(e.c);
        }
        return curr;
    }

    @Override
    public String toString() { 
        String s = "[\n";
        for(Triple<Integer, Integer, Double> t : this.entries){
            s += t.toString()+"\n";
        }
        return s + "]";
    }
    
}
