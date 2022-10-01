package project_utils;

/**
 * Interface for double sequences (e.g., <i>(1,2,3,4,2,3,4,5,3,4,5,6,...)</i>).
 * @author Jonas Schulz
 */
public interface DoubleSequence {
    /**
     * Defines the double of the sequence at the specified index.
     * @param index Parameter for the index in the sequence.
     * @return Value of the sequence at <code>index</code>.
     */
    public double at(int index);
    
    /**
     * Sequence of 1s. (<i>1,1,1,1,1,...</i>)
     */
    public static DoubleSequence one = (int index) -> 1.;
    
    /**
     * Hyperbole sequence <i>a<sub>i</sub> = (i + s) / (i + 1).
     * @param start <i>s</i>.
     * @return The hyperbole sequence with parameter <i>s</i>.
     */
    public static DoubleSequence hyberbole(double start){
        return ((int index) -> (index+start)/(index+1));
    }
    
    /**
     * Constant sequence. (<i>c,c,c,c,c,...</i>)
     * @param c The constant.
     * @return The constant sequence.
     */
    public static DoubleSequence constant(double c){
        return ((int index) -> c);
    }
    
    /**
     * Small example of the use of this class.
     * @param args Program parameters. Currently, none is used.
     */
    public static void main(String[] args){
        DoubleSequence hseq = new DoubleSequence(){
                        public double at(int index){
                            return 1.;
                        }
                    };
        System.out.println("one-constant:");
        System.out.println(hseq.at(0)+", "+hseq.at(1)+", "+hseq.at(2));
        for(double h = 0.4; h <= 1.3; h+= 0.05){
            System.out.println("h is "+h);
            final double h2 = h;
            DoubleSequence hseq2 = new DoubleSequence(){
                        public double at(int index){
                            return h2;
                        }
                    };
            System.out.println(hseq2.at(0)+", "+hseq2.at(1)+", "+hseq2.at(2));
        }
        double h3 = 0.4; for(; h3 < 1; h3 +=0.05);
        System.out.println("h3 is "+h3);
    }
}
