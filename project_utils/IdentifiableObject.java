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
public class IdentifiableObject implements Identifiable{
    Object o;
    final int id;
    private static int static_id = 0;
    public IdentifiableObject(Object o){
        this.o = o;
        this.id = static_id++;
    }

    public Object get(){
        return o;
    }
    @Override
    public int getID() {
        return id;
    }

    @Override
    public String toString() {
        return "#"+id+" "+o.toString();
    }
    
    
}
