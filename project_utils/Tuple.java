package project_utils;

/**
 * A simple class implementing generic binary tuple functionality.
 * @author Jonas Schulz
 * @param <A> Type of the first tuple object.
 * @param <B> Type of the second tuple object.
 */
public class Tuple <A,B> {
    /**
     * First tuple object (a,_).
     */
    public A a;
    /**
     * Second tuple object (_,b).
     */
    public B b;
    /**
     * Standard constructor.
     * @param a First tuple object (a,_).
     * @param b Second tuple object (_,b).
     */
    public Tuple(A a, B b){
        this.a = a;
        this.b = b;
    }
    /**
     * {@inheritDoc }
     */
    @Override
    public String toString() {
        return "("+a.toString()+", "+b.toString()+")";
    }

    /**
     * {@inheritDoc }
     */
    @Override
    public boolean equals(Object obj) {
        if(obj instanceof Tuple){
            Tuple other = (Tuple)obj;
            return other.a.equals(this.a) && other.b.equals(this.b);
        } else {
            return false;
        }
    }
    
}
