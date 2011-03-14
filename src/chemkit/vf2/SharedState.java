package chemkit.vf2;

import java.util.Arrays;

// The SharedState class holds four arrays containing the mapping between
// the two graphs and the terminal sets. It is shared between all the states
// in each isomorphism test.
class SharedState {

    private int[] sourceMapping;
    private int[] targetMapping;
    private int[] sourceTerminalSet;
    private int[] targetTerminalSet;

    public SharedState(int sourceSize, int targetSize) {
        sourceMapping = new int[sourceSize];
        Arrays.fill(sourceMapping, -1);

        targetMapping = new int[targetSize];
        Arrays.fill(targetMapping, -1);

        sourceTerminalSet = new int[sourceSize];
        Arrays.fill(sourceTerminalSet, 0);

        targetTerminalSet = new int[targetSize];
        Arrays.fill(targetTerminalSet, 0);
    }

    @Override
    public String toString() {
        return "src: " + Arrays.toString(getSourceMapping())
                + " trg: " + Arrays.toString(getTargetMapping())
                + " sTS: " + Arrays.toString(getSourceTerminalSet())
                + " tTS: " + Arrays.toString(getTargetTerminalSet());
    }

    /**
     * @return the sourceMapping
     */
    public int[] getSourceMapping() {
        return sourceMapping;
    }
    /**
     * @return the targetMapping
     */
    public int[] getTargetMapping() {
        return targetMapping;
    }

    /**
     * @return the sourceTerminalSet
     */
    public int[] getSourceTerminalSet() {
        return sourceTerminalSet;
    }
    
    /**
     * @return the targetTerminalSet
     */
    public int[] getTargetTerminalSet() {
        return targetTerminalSet;
    }
}
