/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grid_algorithm;

import project_utils.Vector;

/**
 *
 * @author kroka
 */
public interface RoutingGraph {
    public Vector approximateCongestion(Vector flow);
}
