/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package project_utils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author kroka
 */
public class SparseVector extends Vector {
    private int n;
    protected HashMap<Integer, Double> entries = new HashMap<>();
    //private Set<Integer> keys = new HashSet<>();
    //private boolean keys_updated = false;

    public SparseVector(int n){
        this.n = n;
    }

    @Override
    public Iterator<Map.Entry<Integer,Double>> iterator() {
        return entries.entrySet().iterator();
    }
    
    @Override
    public int getN() {
        return n;
    }

    @Override
    public void set(int i, double value) {
        entries.put(i, value);
        //keys_updated = false;
    }
    
    public void add(int i, double value){
        boolean alreadyContained = entries.containsKey(i);
        double sum = (alreadyContained ? entries.get(i) : 0) + value;
        entries.put(i, sum);
        //if(!alreadyContained) keys.add(i);
    }
    public SparseVector add_sparse(SparseVector other){
        //update_keys();
        //other.update_keys();
        SparseVector res = (SparseVector)clone();
        //for(int key : keys){
        //    res.set(key, get(key));
        //}
//        for(int key : other.keys){
//            res.add(key, other.get(key));
//        }
        for(Map.Entry<Integer,Double> e : other){
            res.add(e.getKey(), e.getValue());
        }
        return res;
    }

    @Override
    protected Object clone(){
        SparseVector copy = new SparseVector(n);
        for(Map.Entry<Integer, Double> e : entries.entrySet()){
            copy.set(e.getKey(), e.getValue());
        }
        return copy;
    }
    

    
    @Override
    public double get(int i) {
        return entries.get(i);
    }

//    @Override
//    public Matrix toMatrix() {
//        SparseMatrix sm = new SparseMatrix(getN(), 1);
//        update_keys();
//        for(int key : keys){
//            sm.set(key, 0, get(key));
//        }
//        return sm;
//    }
    

    static boolean sparse_string = false;
    @Override
    public String toString() {
        return sparse_string ? toSparseString() : toCompleteString();
    }
    public String toSparseString(){
        return Arrays.toString(entries.entrySet().toArray());
    }
    public String toCompleteString(){
        return super.toString();
    }
}
