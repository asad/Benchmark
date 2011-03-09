package vf;

/**
 * 
 * @author Asad
 * @param <T>
 * @param <S>
 */
class Pair<T, S> {

    private T first;
    private S second;

    public Pair(T a, S b) {
        this.first = a;
        this.second = b;
    }

    @Override
    public String toString() {
        return "(" + getFirst() + ", " + getSecond() + ")";
    }

    /**
     * @return the first
     */
    public T getFirst() {
        return first;
    }

    /**
     * @param first the first to set
     */
    public void setFirst(T first) {
        this.first = first;
    }

    /**
     * @return the second
     */
    public S getSecond() {
        return second;
    }

    /**
     * @param second the second to set
     */
    public void setSecond(S second) {
        this.second = second;
    }
}