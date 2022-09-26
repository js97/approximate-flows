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
