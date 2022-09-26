package project_utils;

/**
 *
 * @author Jonas Schulz
 */
public interface DoubleSequence {
    public double at(int index);
    
    public static DoubleSequence one = (int index) -> 1.;
    public static DoubleSequence hyberbole(double start){
        return ((int index) -> (index+start)/(index+1));
    }
    public static DoubleSequence constant(double c){
        return ((int index) -> c);
    }
    
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
