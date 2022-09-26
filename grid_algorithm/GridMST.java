package grid_algorithm;

import java.util.Arrays;

/**
 * This class represents a maximal spanning tree for grid graphs.
 * The construction is described within the master's thesis.
 * Currently, it only supports the routing algorithm.
 * @author Jonas Schulz
 */
public class GridMST {
    /**
     * The grid graph related to this maximum spanning tree.
     */
    GridGraph g;
    /**
     * Standard constructor.
     * @param g The grid graph on which the maximum spanning tree should be constructed.
     */
    public GridMST(GridGraph g){
        this.g = g;
    }
    /**
     * Calculates the next vertex in the traversal, based on the iteration scheme.
     * The iteration scheme follows a valid leaf-elimination algorithm, s.t.
     * this method only returns vertices whose predecessors have been eliminated,
     * when traversed from the start vertex. <br>
     * Example: with graph dimensions (2,5,3), the start index should be [1,4,2].
     * Then, we would get: 
     * <ul>
     * <li><code>iter_next([1,4,2])<code> = <code>[1,4,1]</code>,</li>
     * <li><code>iter_next([1,4,1])<code> = <code>[1,4,0]</code>,</li>
     * <li><code>iter_next([1,4,0])<code> = <code>[1,3,2]</code>,</li>
     * <li><code>iter_next([1,3,2])<code> = <code>[1,3,1]</code>,</li>
     * <li>...,</li>
     * <li><code>iter_next([0,0,1])<code> = <code>[0,0,0]</code>,</li>
     * <li><code>iter_next([0,0,0])<code> = <code>null</code>.</li>
     * </ul>
     * @param current The current index in vector grid coordinates.
     * @return The vector grid coordinates of the vertex to be eliminated after <code>current</code>.
     */
    public int[] iter_next(int[] current){
        return iter_next(g, current);
    }
    /**
     * Implements the routing of a demand through this MST.
     * For more details of the algorithm, look at the master's thesis. For more
     * details regarding the implementation, also see the comments in the code.
     * @param b Demand that should be routed through the MST of <code>G</code>.
     * @return The flow on this MST that satisfies demand <code>b</code>.
     */
    public GridFlow route(GridDemand b){
        return route(g, b);
    }
    /**
     * Optimized version of <code>route</code>.
     * Instead of iterating the dimensions and using vector grid coordinates,
     * this algorithm uses the enumeration scheme of the vertices. The vertex traversal
     * remains the same.
     * @param b Demand that should be routed through this MST.
     * @return The flow on this MST that satisfies demand <code>b</code>.
     */
    public GridFlow route2(GridDemand b){
        return route2(g, b);
    }
    /**
     * Implements the routing of a demand through a maximal spanning tree of a grid graph.
     * For more details of the algorithm, look at the master's thesis. For more
     * details regarding the implementation, also see the comments in the code.
     * @param G Grid graph on which the MST should be constructed.
     * @param b Demand that should be routed through the MST of <code>G</code>.
     * @return The flow on the MST of <code>G</code> that satisfies demand <code>b</code>.
     */
    public static GridFlow route(GridGraph G, GridDemand b){
        GridFlow f = new GridFlow(G);
        GridDemand r = b.clone();
        int[] n = G.nodesPerDim;
        int d = n.length;
        
        /**
         * We iterate the dimensions. At iteration i, we handle all edges of the
         * MST that are directed along this dimension, e.g. the edge [1,1,2,2,1] -- [1,2,2,2,1]
         * would be handled in iteration i = 1. For a more detailed iteration scheme,
         * look at the master's thesis.
         */
        for (int i = 0; i < d; i++) {
            /**
             * If n[i] == 1, there is no iteration in this dimension.
             */
            if(n[i] == 1) continue;
            /**
             * Set start and end vertices.
             */
            int[] v_start = new int[d];
            int[] v_end = new int[d];
            for (int j = 0; j < d; j++) {
                if(j < i){
                    v_start[j] = 0;
                    v_start[j] = 0;
                } else if(j == i){
                    v_start[j] = n[j] - 1;
                    v_end[j] = 1;
                } else {
                    v_start[j] = n[j] - 1;
                    v_end[j] = 0;
                }
            }
            /**
             * Iterate from start to end with step iter_next.
             * The end is reached iff the arrays v and v_end equal.
             */
            boolean last_reached = false;
            for(int[] v = v_start; !last_reached; v = iter_next(G, v)){
                
                int[] u = Arrays.copyOf(v, v.length);
                /**
                 * The maximal spanning tree we construct will route from this vertex (v)
                 * to u. v is a leaf of the remaining subgraph, it is eliminated here.
                 * The neighbour in this dimension i is u. In accordance to above, we just
                 * set u = v, but u[i] = v[i] - 1, and update the residual demand.
                 */
                u[i] = v[i]-1;
                f.set(u, v, r.get(v));
                r.set(u, r.get(v)+r.get(u));
                r.set(v, 0.);
                if(Arrays.equals(v, v_end)){
                    last_reached = true;
                }
            }
            
        }
        
        return f;
    }
    /**
     * Calculates the next vertex in the traversal, based on the iteration scheme.
     * The iteration scheme follows a valid leaf-elimination algorithm, s.t.
     * this method only returns vertices whose predecessors have been eliminated,
     * when traversed from the start vertex. <br>
     * Example: with <code>g = new GridGraph(2,5,3)</code>, the start index should be [1,4,2].
     * Then, we would get: 
     * <ul>
     * <li><code>iter_next(g, [1,4,2])<code> = <code>[1,4,1]</code>,</li>
     * <li><code>iter_next(g, [1,4,1])<code> = <code>[1,4,0]</code>,</li>
     * <li><code>iter_next(g, [1,4,0])<code> = <code>[1,3,2]</code>,</li>
     * <li><code>iter_next(g, [1,3,2])<code> = <code>[1,3,1]</code>,</li>
     * <li>...</li>,
     * <li><code>iter_next(g, [0,0,1])<code> = <code>[0,0,0]</code>,</li>
     * <li><code>iter_next(g, [0,0,0])<code> = <code>null</code>.</li>
     * </ul>
     * @param G The grid graph on which the MST should be constructed on.
     * @param current The current index in vector grid coordinates.
     * @return The vector grid coordinates of the vertex to be eliminated after <code>current</code>.
     */
    static int[] iter_next(GridGraph G, int[] current){
        int last_nonzero = current.length-1;
        while(last_nonzero >= 0 && current[last_nonzero] == 0){
            last_nonzero--;
        }
        int[] next = Arrays.copyOf(current, current.length);
        if(last_nonzero >= 0){
            next[last_nonzero]--;
            for (int i = last_nonzero+1; i < current.length; i++) {
                next[i] = G.nodesPerDim[i]-1;
            }
            return next;
        } else {
            return null;
        }
    }
    /**
     * Optimized version of <code>route</code>.
     * Instead of iterating the dimensions and using vector grid coordinates,
     * this algorithm uses the enumeration scheme of the vertices. The vertex traversal
     * remains the same.
     * @param G Grid graph on which the MST should be constructed.
     * @param b Demand that should be routed through the MST of <code>G</code>.
     * @return The flow on the MST of <code>G</code> that satisfies demand <code>b</code>.
     */
    public static GridFlow route2(GridGraph G, GridDemand b){
        GridFlow f = new GridFlow(G);
        GridDemand r = b.clone();
        int n = G.getN();
        int i = 0;
        for (int k = n - 1; k >= 1; k--) {
            int[] v = G.toPosition(k);
            /**
             * The edge of the currently eliminated leaf v can be determined as
             * in route(), which is why we calculate it with vertex grid coordinates.
             * It could be done more efficiently with subtracting n_(i+1)*...*n_d, but
             * this would be worse code design, as it needs to be consistent with the implementation
             * of the toIndex function of the GridGraph class.
             */
            while(v[i] == 0){
                i++;
            }
            v[i]--;
            int u = G.toIndex(v);
            f.set(u, k, r.get(k));
            r.set(u, r.get(u)+r.get(k));
            r.set(k, 0.);
        }
        return f;
    }
}
