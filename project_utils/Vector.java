/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package project_utils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.lang.Math;
import java.util.AbstractMap;
import java.util.Iterator;
import java.util.Map.Entry;

/**
 *
 * @author kroka
 */
public abstract class Vector implements Iterable<Entry<Integer,Double>>{
    public abstract int getN();
    public abstract void set(int i, double value);
    public abstract double get(int i);
    public double scalar(Vector other){
        double sum = 0.;
        for (Entry<Integer, Double> e : this) {
            sum += e.getValue()*other.get(e.getKey());
        }
        return sum;
    }

    //TODO: change "get(e.getKey())" in loops over thie (Entry e) to
    // while(this.iterator.hasNext() && other.iterator.hasNext())
    // and simultanously call next() on both iterators when calculation
    // component-wise to optimize runtime
    // (in addition, only perform action if e1.key==e2.key and if not,
    // only call e[min(e1.key,e2.key)].next())
    
    @Override
    public Iterator<Entry<Integer,Double>> iterator() {
        return new Iterator<Entry<Integer,Double>>() {
            int i = 0;
            @Override
            public boolean hasNext() {
                return i < getN();
            }

            @Override
            public Entry<Integer,Double> next() {
                Double next = get(i);
                i++;
                return new AbstractMap.SimpleEntry<>(i-1,next);
            }
        };
    }
    
//    public Matrix toMatrix() {
//        double[][] res = new double[1][getN()];
//        for(int i = 0; i < getN(); i++){
//            res[0][1] = get(i);
//        }
//        return new ArrayMatrix(res);
//    }
    public Vector scalarMultiply(double s){
        double[] res = new double[getN()];
        for(Entry<Integer,Double> e : this){
            res[e.getKey()] = s*e.getValue();
        }
        return new ArrayVector(res);
    }
    public double Linf(){
        double max = get(0);
        for(Entry<Integer,Double> e : this){
            double val = Math.abs(e.getValue());
            if(val > max) max = val;
        }
        return max;
    }
    public double L1(){
        double sum = 0.;
        for(Entry<Integer,Double> e : this){
            sum += Math.abs(e.getValue());
        }
        return sum;
    }
    public double L2_sq(){
        double sqsum = 0.;
        for(Entry<Integer,Double> e : this){
            sqsum += e.getValue()*e.getValue();
        }
        return sqsum;
    }
    public double L2(){
        return Math.sqrt(L2_sq());
    }
    public Vector add(Vector other){
        Vector res = new ArrayVector(getN());
        for(Entry<Integer,Double> e : this){
            res.set(e.getKey(), e.getValue()+other.get(e.getKey()));
        }
        return res;
    }
    public double softmax(){
        return Math.log(expsmax());
    }
    public double expsmax(){
        double sum = 0;
        int non_zeros = 0;
        for(Entry<Integer,Double> e : this){
            sum += Math.exp(e.getValue());
            non_zeros++;
        }
        sum += getN() - non_zeros;
        return sum;
    }
    public double symm_softmax(){
        return Math.log(explmax());
    }
    public double explmax(){
        double sum = 0;
        int non_zeros = 0;
        for(Entry<Integer,Double> e : this){
            double entr = e.getValue();
            sum += Math.exp(entr);
            sum += Math.exp(-1*entr);
            non_zeros++;
        }
        sum += 2 * (getN() - non_zeros);
        return sum;
    }
    public Vector grad_symm_softmax(){
        double divisor = explmax();
        double[] entr = new double[getN()];
        for(int key = 0; key < getN(); key++){
            entr[key] = Math.exp(get(key))-Math.exp((-1)*get(key));
            entr[key] /= divisor;
        }
        return new ArrayVector(entr);
    }
    @Override
    public String toString() {
        return Arrays.toString(toArray());
    }
    public double[] toArray(){
        double[] entries = new double[getN()];
        for(int i = 0; i < getN(); i++){ 
            entries[i] = get(i);
        }
        return entries;
    }
}
