package vf2;

import java.util.Arrays;

// The SharedState class holds four arrays containing the mapping between
// the two graphs and the terminal sets. It is shared between all the states
// in each isomorphism test.
class SharedState {

    int[] sourceMapping;
    int[] targetMapping;
    int[] sourceTerminalSet;
    int[] targetTerminalSet;

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
        return "src: " + Arrays.toString(sourceMapping)
                + " trg: " + Arrays.toString(targetMapping)
                + " sTS: " + Arrays.toString(sourceTerminalSet)
                + " tTS: " + Arrays.toString(targetTerminalSet);
    }
}
