package project_utils;

import java.util.Arrays;

/**
 *
 * @author kroka
 */
public class ArrayVector extends Vector{
    double[] entries;
    int n;
    
    public ArrayVector(int n){
        entries = new double[n];
        this.n = n;
    }
    public ArrayVector(double[] entries){
        this.entries = entries;
        this.n = entries.length;
    }
    
    public void set(int i, double value){
        entries[i] = value;
    }
    public double get(int i){
        return entries[i];
    }

    @Override
    public int getN() {
        return entries.length;
    }

    @Override
    public String toString() {
        return Arrays.toString(entries);
    }

    @Override
    protected Object clone() {
        return new ArrayVector(super.toArray());
    }
    
    
    
}
