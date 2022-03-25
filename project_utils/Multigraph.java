/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package project_utils;

import java.util.ArrayList;
import java.util.Optional;

/**
 *
 * @author kroka
 */
public class Multigraph {
    ArrayList<ArrayList<IdentifiableObject>> multiVertices;
    ArrayList<ArrayList<Tuple<Integer, Double>>> edges;
    public Multigraph(){
        multiVertices = new ArrayList<>();
        edges = new ArrayList<>();
    }
    
    public void addMultiVertex(){
        multiVertices.add(new ArrayList<>());
    }
    public void addMultiVertex(IdentifiableObject... contained){
        ArrayList<IdentifiableObject> a = new ArrayList<>();
        for(IdentifiableObject io : contained){
            a.add(io);
        }
        multiVertices.add(a);
    }
    public void contract(int multiVertexIndexInto, int multiVertexIndexFrom){
        
    }
    
    public Optional<Integer> getIndexOfObject(IdentifiableObject o){
        for(int i = 0; i < multiVertices.size(); i++){
            for(IdentifiableObject io : multiVertices.get(i)){
                if(Identifiable.equals(o, io)){
                    return Optional.of(i);
                }
            }
        }
        return Optional.empty();
    }
}
