package vf2;

import org.openscience.cdk.annotations.TestClass;

/**
 * Interface for the storing the states of the mapping in the VF algorithm.
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.algorithm.vflib.VFLibTest")
public interface IState {

    /**
     * Returns the current mapping of query atoms onto target atoms.
     * This map is shared among all states obtained through nextState.
     *
     * @return the current mapping of query atoms onto target atoms
     */
    public AtomMapping getMapping();

    /**
     * Returns true if another candidate match can be found or
     * false otherwise.
     *
     * @param candidate 
     * @return true if another candidate mapping can be found or
     * false otherwise.
     */
    public boolean hasNextCandidate(Match<Integer, Integer> candidate);

    /**
     * Returns the next candidate match.
     *
     * @param lastCandidate 
     * @return the next candidate match.
     */
    public Match<Integer, Integer> nextCandidate(Match<Integer, Integer> lastCandidate);

    /**
     * Returns true if the given match will work with the current
     * map, or false otherwise.
     *
     * @param match the match to consider
     * @return true if the given match will work with the current
     * map, or false otherwise.
     */
    public boolean isMatchFeasible(Match<Integer, Integer> match);

    /**
     * Returns true if all atoms in the query molecule have been
     * mapped.
     *
     * @return true if all atoms in the query molecule have been
     * mapped.
     */
    public boolean isGoal();

    /**
     * Returns true if no match will come from this IState.
     *
     * @return true if no match will come from this IState
     */
    public boolean isDead();

    /**
     * Potential candidates are added
     * to the current mapping.
     *
     * @param candidate 
     */
    public void addPair(Match<Integer, Integer> candidate);

    /**
     * Returns this IState's atom map to its original condition.
     */
    public void backTrack();
}
