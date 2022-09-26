package project_utils;

/**
 *
 * @author kroka
 */
public class Tuple <A,B> {
    public A a;
    public B b;
    public Tuple(A a, B b){
        this.a = a;
        this.b = b;
    }
    @Override
    public String toString() {
        return "("+a.toString()+", "+b.toString()+")";
    }

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
