/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package project_utils;

/**
 *
 * @author kroka
 */
public class Triple<A, B, C> {
    public A a;
    public B b;
    public C c;
    public Triple(A a, B b, C c){
        this.a = a;
        this.b = b;
        this.c = c;
    }

    @Override
    public String toString() {
        return "("+a.toString()+", "+b.toString()+", "+c.toString()+")";
    }
    @Override
    public boolean equals(Object obj) {
        if(obj instanceof Triple){
            Triple other = (Triple)obj;
            return other.a.equals(this.a) && other.b.equals(this.b) && other.c.equals(this.c);
        } else {
            return false;
        }
    }

    
}
