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
public class Math {
    public static Tuple<Integer, Double> max(double... vals){
        double max = Double.MIN_VALUE;
        int index = -1;
        int cur = 0;
        for(double d : vals){
            if(d > max){
                max = d;
                index = cur;
            }
            cur++;
        }
        return new Tuple<>(index, max);
    }
}
