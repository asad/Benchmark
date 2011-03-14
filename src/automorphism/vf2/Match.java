package automorphism.vf2;

/**
 * 
 * @author Asad
 * @param <T>
 * @param <S>
 */
public class Match<T, S> {

    private T source;
    private S target;

    public Match(T a, S b) {
        this.source = a;
        this.target = b;
    }

    @Override
    public String toString() {
        return "(" + getSourceAtom() + ", " + getTargetAtom() + ")";
    }

    /**
     * @return the source
     */
    public T getSourceAtom() {
        return source;
    }

    /**
     * @param first 
     * @param source the source to set
     */
    public void setSourceAtom(T first) {
        this.source = first;
    }

    /**
     * @return the target
     */
    public S getTargetAtom() {
        return target;
    }

    /**
     * @param second 
     * @param target the target to set
     */
    public void setTargetAtom(S second) {
        this.target = second;
    }
}