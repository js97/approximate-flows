/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grid_algorithm;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map.Entry;
import project_utils.Triple;
import project_utils.Tuple;

/**
 *
 * @author kroka
 */
public class HashedGridFlow {
    GridGraph g;
    HashMap<Integer, Double> entries;
    int m;
    
    public HashedGridFlow(GridGraph g){
        this.g = g;
        calculateM();
        entries = new HashMap<>(m);
    }
    
    HashedGridFlow(HashedGridFlow f, boolean clone){
        this.entries = new HashMap<>();
        this.g = f.g;
        this.m = f.m;
        if(clone){
            for(Entry<Integer, Double> t : f.entries.entrySet()){
                this.entries.put(t.getKey(), t.getValue());
            }
        }
    }
    
    public int calculateM(){
        int sum = 0;
        for(int i = 0; i < g.nodesPerDim.length; i++){
            sum += (g.getN()/g.nodesPerDim[i])*(g.nodesPerDim[i]-1);
        }
        // 2*sum or sum ? TODO
        this.m = sum;
        return this.m;
    }

    public int getM() {
        return m;
    }

    public int getN() {
        return g.getN();
    }
    
    public void set(Integer from, Integer to, Double value){
        Triple<Integer, Integer, Double> e = (from > to) ? new Triple<>(to, from, -value) : new Triple<>(from, to, value);
        entries.put(this.hashEntry(new Tuple<>(e.a, e.b)), e.c);
    }
    public void set(int[] from, int[] to, Double value){
        int from_int = g.toIndex(from);
        int to_int = g.toIndex(to);
        this.set(from_int, to_int, value);
    }
    public void add(Integer from, Integer to, double value){
        if(from > to)
            this.add(this.hashEntry(new Tuple<>(to, from)), -value);
        else
            this.add(this.hashEntry(new Tuple<>(from, to)), value);
    }
    void set(int combinedIndex, double value){
        entries.put(combinedIndex, value);
    }
    void add(int combinedIndex, double value){
        entries.put(combinedIndex, (entries.containsKey(combinedIndex) ? entries.get(combinedIndex) : 0)+ value);
    }
    
    public double get(Integer from, Integer to){
        int combined = this.hashEntry((from > to) ? (new Tuple<>(to, from)) : (new Tuple<>(from, to)));
        return (entries.containsKey(combined)) ? ((from > to ? -1 : 1)*entries.get(combined)) : 0.;
    }
    public double get(int[] from, int[] to){
        return get(this.g.toIndex(from), this.g.toIndex(to));
    }
    
    // this method represents the calculation B*f for 
    // - B the incidence matrix
    // - f the flow
    public GridDemand calculateExcessFlows(){
        GridDemand b = new GridDemand(this);
        for(Entry<Integer, Double> e : entries.entrySet()){
            Tuple<Integer, Integer> edge = this.hashToTuple(e.getKey());
            b.add(edge.a, -e.getValue());
            b.add(edge.b, e.getValue());
        }
        return b;
    }
    
    public double lmax_exp(){
        double sum = 0;
        int non_zeros = 0;
        for(Entry<Integer, Double> e : entries.entrySet()){
            double exp = Math.exp(e.getValue());
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
    public double lmax_exp_shifted(double shift){
        double sum_exps = 0;
        int nonzeros = 0;
//        System.out.println("lmax_exp_shifted: ");
        for(Entry<Integer, Double> e : entries.entrySet()){
            double exp = Math.exp(e.getValue() + shift) + Math.exp((-e.getValue()) + shift);
            sum_exps += exp;
            nonzeros++;
//            System.out.println("    edge value: "+e.getValue());
        }
//        System.out.println("  shift is "+shift);
//        System.out.println("  sum_exps is "+sum_exps);
        double exp_shift = Math.exp(shift);
//        System.out.println("  exp_shift is "+exp_shift);
        sum_exps += 2 * exp_shift * (m - nonzeros);
//        System.out.println("  there are "+nonzeros+" nonzeros and "+(m-nonzeros)+" zeros.");
        return sum_exps;
    }
    public double lmax(){
        // might use shifted version with shift = 0
        return Math.log(lmax_exp());
    }
    public double lmax_shifted(double shift){
        return Math.log(lmax_exp_shifted(shift)) - shift;
    }
    public double lmax_shifted(){
        double shift = getDefaultShift();
        return lmax_shifted(shift);
    }
    public HashedGridFlow gradient_lmax(){
        double sum_exps = lmax_exp();
        HashedGridFlow f = new HashedGridFlow(this, false);
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            double exp = Math.exp(t.getValue());
            f.set(t.getKey(), (exp-(1./exp))/sum_exps);
        }
        // where there is i with no entry, xi = 0, hence e^xi - e^(-xi) = 0, thus 
        // also the gradient is 0 at i
        return f;
    }
    public HashedGridFlow gradient_lmax_shifted(double shift){
//        System.out.println("Calculating gradient with shift "+shift);
        double sum_exps = lmax_exp_shifted(shift);
//        System.out.println("sum_exps is "+sum_exps);
//        System.out.println("Maximal entry is "+this.linf());
        HashedGridFlow f = new HashedGridFlow(this, false);
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            double exp1 = Math.exp(t.getValue() + shift);
            double exp2 = Math.exp(-t.getValue() + shift);
            f.set(t.getKey(), (exp1 - exp2) / sum_exps);
//            System.out.println("exponent: "+t.getValue());
//            System.out.println("exp1: "+exp1);
//            System.out.println("exp2: "+exp2);
        }
//        System.out.println("Gradient l1 is "+f.l1()+" and linf is "+f.linf());
        return f;
    }
    double getDefaultShift(){
        return -1*linf();
    }
    public HashedGridFlow gradient_lmax_shifted(){
        return gradient_lmax_shifted(getDefaultShift());
    }
    public double linf(){
        double max = 0;
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            double v = Math.abs(t.getValue());
            if(v > max) max = v;
        }
        return max;
    }

    public HashedGridFlow scale(double d) {
        HashedGridFlow cl = new HashedGridFlow(this, true);
        for(Entry<Integer, Double> t : cl.entries.entrySet()){
            t.setValue(t.getValue() * d);
        }
        return cl;
    }
    public HashedGridFlow scale_inplace(double d) {
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            t.setValue(t.getValue() * d);
        }
        return this;
    }
    public static HashedGridFlow add(HashedGridFlow a, HashedGridFlow b) {
        HashedGridFlow cl = new HashedGridFlow(a, true);
//        for(Entry<Integer, Double> t : a.entries.entrySet()){
//            cl.set(t.getKey(), t.getValue());
//        }
        for(Entry<Integer, Double> t : b.entries.entrySet()){
            cl.add(t.getKey(), t.getValue());
        }
        // cl.contractDuplicates(true);
        return cl;
    }
    
    int hashEntry(Tuple<Integer, Integer> t){
        return (((int)t.a) << 16) + t.b;
    }
    Tuple<Integer, Integer> hashToTuple(int hash){
        int second = hash & 0x0000_ffff;
        int first = (hash >> 16) & 0x0000_ffff;
        return new Tuple<>(first, second);
    }
    
    public double l1() {
        double curr = 0;
        for(Entry<Integer, Double> e : entries.entrySet()){
            curr += Math.abs(e.getValue());
        }
        return curr;
    }

    @Override
    public String toString() { 
        String s = "[\n";
        for(Entry<Integer, Double> t : this.entries.entrySet()){
            Tuple<Integer, Integer> i = this.hashToTuple(t.getKey());
            s += "\t"+i.toString()+" :   "+t.getValue()+"\n";
        }
        return s + "]";
    }
    
    public String tikz2D(){
//        String s = "";
        String s = "\\begin{tikzpicture}[roundnode/.style={circle, draw=green!60, fill=green!5, very thick, minimum size=7mm}, scale=1.2]\n";
        //String s = "\\node (anchor) {};\n";
        DecimalFormat df = new DecimalFormat("#.000");
        for(int i = 0; i < getN(); i++){
            double scale = 2;
            s += "\\node[roundnode, minimum size = 1cm] at ("+scale*g.toPosition(i)[0]+", "+scale*g.toPosition(i)[1]+") ("+i+") {\\textcolor{black}{"+i+"}};\n";
        }
        for (int i = 0; i < getN(); i++) {
            int[][] neighbours = g.getNeighbours(g.toPosition(i));
            for(int[] n : neighbours){
                int indexTo = g.toIndex(n);
                if(indexTo > i){
                    String number = (this.get(i, indexTo) == 0.0) ? "\\textcolor{blue!28}{0}" : df.format(this.get(i, indexTo));
                    if(indexTo == i + 1)
                        s += "\\draw[thick, ->] ("+i+") -- ("+indexTo+") node[midway,right] {"+number+"};\n";
                    else
                        s += "\\draw[thick, ->] ("+i+") -- ("+indexTo+") node[midway,above] {"+number+"};\n";
                }
            }
            
        }
        s += "\\end{tikzpicture}";
        return s;
    }
}
