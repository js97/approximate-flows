package project_utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import jdk.jshell.spi.ExecutionControl;

/**
 *
 * @author kroka
 */
public class Graph {
    Map<Integer, IdentifiableObject> node_to_object_map;
    
    //might just use size since V is equal to [V.size()-1]_0
    protected ArrayList<Integer> V;
    protected ArrayList<ArrayList<Tuple<Integer, Double>>> E;
    public Graph(){
        V = new ArrayList<>();
        E = new ArrayList<>();
        node_to_object_map = new HashMap<>();
    }
    // internal addition only, object map and object have to be added afterwards
    private int addNode(){
        int new_node = V.size();
        V.add(new_node);
        E.add(new ArrayList<>());
        return new_node;
    }
    public int addNode(IdentifiableObject o){
        int new_node = addNode();
        node_to_object_map.put(new_node, o);
        return new_node;
    }
    public int addNode(Object o){
        IdentifiableObject obj = new IdentifiableObject(o);
        int new_node = addNode(obj);
        return new_node;
    }
    public void addEdge(int from, int to, double weight){
        if(from < 0 || to < 0 || from >= V.size() || to >= V.size()){
            throw new IllegalArgumentException("Tried to add an edge between non-existent nodes.");
        }
        E.get(from).add(new Tuple<Integer, Double>(to, weight));
    }
    public void addOrUpdateEdge(int from, int to, double weight){
        if (edgeExists(from, to))
            updateEdge(from, to, weight);
        else
            addEdge(from, to, weight);
    }
    public boolean edgeExists(int from, int to){
        int index = getIndexOfEdgeInAdjacencyList(from, to);
        return (index < E.get(from).size());
    }
    protected int getIndexOfEdgeInAdjacencyList(int from, int to){
        ArrayList<Tuple<Integer, Double>> adjList = E.get(from);
        int index = 0;
        for(Tuple<Integer, Double> t : adjList){
            if(t.a.equals(to)){
                return index;
            }
            index++;
        }
        return index;
    }
    public void updateEdge(int from, int to, double new_weight){
        if(from < 0 || to < 0 || from >= V.size() || to >= V.size()){
            throw new IllegalArgumentException("Tried to update an edge between non-existent nodes.");
        }
        if(!edgeExists(from,to)){
            throw new IllegalArgumentException("Tried to update a non-existent edge.");
        }
        int index = getIndexOfEdgeInAdjacencyList(from, to);
        E.get(from).get(index).b = new_weight;
    }
    public Optional<Double> getEdge(int from, int to){
        if(edgeExists(from, to)){
            return Optional.of(E.get(from).get(getIndexOfEdgeInAdjacencyList(from, to)).b);
        } else {
            return Optional.empty();
        }
    }
    public int n(){
        return V.size();
    }
    public int m(){
        // could do edge counter for constant time
        int size = 0;
        for(ArrayList<Tuple<Integer, Double>> adjList : E){
            size += adjList.size();
        }
        return size;
    }

    @Override
    public String toString() {
        String s = "Nodes: {";
        for(Integer i : V){
            if((int)i == 0){
                s += i;
            } else {
                s += ", "+i;
            }
        }
        s += "}\n";
        s += "Edges (Adjacency Lists): {";
        int index = 0;
        for(ArrayList<Tuple<Integer, Double>> adjList : E){
            if(index == 0){
                s += "\n";
            } 
            s += "\t" + index + ":";
            if(adjList.isEmpty()){
                s += "\t" + "no edges";
            } else {
                for(Tuple<Integer, Double> t : adjList){
                    s += "\n\t\t";
                    s += t.toString();
                }
            }
            s += "\n";
            index++;
        }
        s += "}";
        return s;
    }

    @Override
    protected Object clone() {
        ArrayList<Integer> nodes = new ArrayList<>();
        for(int i : V) nodes.add(i);
        ArrayList<ArrayList<Tuple<Integer, Double>>> edges = new ArrayList<>();
        for(ArrayList<Tuple<Integer, Double>> adjList : E){
            ArrayList<Tuple<Integer, Double>> newList = new ArrayList<>();
            for(Tuple<Integer, Double> t : adjList){
                newList.add(new Tuple<>((int) t.a, (double) t.b));
            }
            edges.add(newList);
        }
        Graph cl = new Graph();
        cl.V = nodes;
        cl.E = edges;
        return cl;
    }
    
    public Graph contract_ext(int a, int b) {
        Graph g = new Graph();
        for (int i = 0; i < V.size() - 1; i++) {
            g.addNode();
            // new graph doesn't contain b
            // index b will be taken from last element
            if(i == b){
                g.node_to_object_map.put(i, this.node_to_object_map.get(V.size()-1));
            } else {
                g.node_to_object_map.put(i, this.node_to_object_map.get(i));
            }
        }
        // if b == V.size() - 1, nothing has been and also has to be done for vertices
        int cfrom = 0;
        for (ArrayList<Tuple<Integer, Double>> adjList : E){
            for (Tuple<Integer, Double> e : adjList){
                // edges from/to b have to be merged with edges from/to a
                // edge between a and b is to be removed
                // all other edges can be copied
                // -> new edge has to use index a instead of index b (if merge)
                // -> new edge has to use index b instead of index last
                int cto = (int)e.a;
                int nfrom = cfrom;
                if(cfrom == b){
                    nfrom = a;
                } else if(cfrom == V.size()-1){
                    nfrom = b;
                }
                int nto = cto;
                if(cto == b){
                    nto = a;
                } else if(cto == V.size() - 1){
                    nto = b;
                }
                if(nfrom != nto){
                    g.addEdge(nfrom, nto, e.b);
                }
            }
            cfrom++;
        }
        return g;
    }
    public Optional<Integer> translate_id_from(Graph from, int id){
        IdentifiableObject obj = from.node_to_object_map.get(id);
        for(Integer v : V){
            if(Identifiable.equals(node_to_object_map.get(v), obj)){
                return Optional.of(v);
            }
        }
        return Optional.empty();
    }
}
