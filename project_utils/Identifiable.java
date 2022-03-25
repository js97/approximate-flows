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
public interface Identifiable {
    public int getID();
    public static boolean equals(Identifiable a, Identifiable b){
        return a.getID() == b.getID();
    }
}
