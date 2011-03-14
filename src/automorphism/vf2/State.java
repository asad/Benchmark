package automorphism.vf2;

import automorphism.vf2.matcher.VFBondMatcher;
import automorphism.vf2.matcher.AtomMatcher;
import automorphism.vf2.matcher.VFAtomMatcher;
import automorphism.vf2.matcher.BondMatcher;
import java.util.List;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

// The State class represents a single state in the isomorphism detection
// algorithm. Every state uses and modifies the same SharedState object.
class State implements IState {

    int getSize() {
        return size;
    }

    IAtomContainer getSource() {
        return source;
    }

    IAtomContainer getTarget() {
        return target;
    }

    IAtom sourceAtom(int index) {
        return source.getAtom(index);
    }

    IAtom targetAtom(int index) {
        return target.getAtom(index);
    }
    int size;
    int sourceTerminalSize;
    int targetTerminalSize;
    IAtomContainer source;
    IAtomContainer target;
    Match<Integer, Integer> candidates;
    SharedState sharedState;
    boolean ownSharedState;

    State(IAtomContainer source, IAtomContainer target) {
        this.size = 0;
        this.sourceTerminalSize = 0;
        this.targetTerminalSize = 0;
        this.source = source;
        this.target = target;
        this.candidates = new Match<Integer, Integer>(-1, -1);
        this.sharedState = new SharedState(source.getAtomCount(),
                target.getAtomCount());
        this.ownSharedState = true;
    }

    State(State state) {
        this.size = state.size;
        this.sourceTerminalSize = state.sourceTerminalSize;
        this.targetTerminalSize = state.targetTerminalSize;
        this.source = state.source;
        this.target = state.target;
        this.candidates = new Match<Integer, Integer>(-1, -1);
        this.sharedState = state.sharedState;
        this.ownSharedState = false;
    }

    // Returns true if the state contains an isomorphism.
    @Override
    public boolean isGoal() {
        return size == source.getAtomCount();
    }

    @Override
    public boolean isDead() {
        return source.getAtomCount() > target.getAtomCount();
    }

    @Override
    public boolean hasNextCandidate(Match<Integer, Integer> candidate) {
        return candidate.getSourceAtom() == -1 ? false : true;
    }

    // Returns the current isomorphism for the state in an AtomMapping
    // object.
    @Override
    public AtomMapping getMapping() {
        AtomMapping mapping = new AtomMapping(source, target);

        for (int i = 0; i < size; i++) {
            mapping.add(source.getAtom(i),
                    target.getAtom(sharedState.getSourceMapping()[i]));
        }
        return mapping;
    }

    // Returns the next candidate pair (sourceAtom, targetAtom) to be added
    // to the state. The candidate should be checked for feasibility and then added
    // using the addPair() method.
    @Override
    public Match<Integer, Integer> nextCandidate(
            Match<Integer, Integer> lastCandidate) {
        int lastSourceAtom = lastCandidate.getSourceAtom();
        int lastTargetAtom = lastCandidate.getTargetAtom();

        int sourceSize = source.getAtomCount();
        int targetSize = target.getAtomCount();

        if (lastSourceAtom == -1) {
            lastSourceAtom = 0;
        }

        if (lastTargetAtom == -1) {
            lastTargetAtom = 0;
        } else {
            lastTargetAtom++;
        }

        if (sourceTerminalSize > size && targetTerminalSize > size) {
            while (lastSourceAtom < sourceSize
                    && (sharedState.getSourceMapping()[lastSourceAtom] != -1 || sharedState.getSourceTerminalSet()[lastSourceAtom] == 0)) {
                lastSourceAtom++;
                lastTargetAtom = 0;
            }
        } else {
            while (lastSourceAtom < sourceSize
                    && sharedState.getSourceMapping()[lastSourceAtom] != -1) {
                lastSourceAtom++;
                lastTargetAtom = 0;
            }
        }

        if (sourceTerminalSize > size && targetTerminalSize > size) {
            while (lastTargetAtom < targetSize
                    && (sharedState.getTargetMapping()[lastTargetAtom] != -1 || sharedState.getTargetTerminalSet()[lastTargetAtom] == 0)) {
                lastTargetAtom++;
            }
        } else {
            while (lastTargetAtom < targetSize
                    && sharedState.getTargetMapping()[lastTargetAtom] != -1) {
                lastTargetAtom++;
            }
        }

        if (lastSourceAtom < sourceSize && lastTargetAtom < targetSize) {
            return new Match<Integer, Integer>(lastSourceAtom, lastTargetAtom);
        }

        return new Match<Integer, Integer>(-1, -1);
    }

    // Adds the candidate pair (sourceAtom, targetAtom) to the state. The
    // candidate pair must be feasible to add it to the state.
    @Override
    public void addPair(Match<Integer, Integer> candidate) {
        size++;
        candidates = candidate;

        int sourceAtom = candidate.getSourceAtom();
        int targetAtom = candidate.getTargetAtom();

        if (sharedState.getSourceTerminalSet()[sourceAtom] < 1) {
            sharedState.getSourceTerminalSet()[sourceAtom] = size;
//                sourceTerminalSize++;
        }

        if (sharedState.getTargetTerminalSet()[targetAtom] < 1) {
            sharedState.getTargetTerminalSet()[targetAtom] = size;
//                targetTerminalSize++;
        }

        sharedState.getSourceMapping()[sourceAtom] = targetAtom;
        sharedState.getTargetMapping()[targetAtom] = sourceAtom;

        List<IAtom> sourceNeighbours =
                source.getConnectedAtomsList(source.getAtom(sourceAtom));
        for (IAtom neighbor : sourceNeighbours) {
            int neighbourIndex = source.getAtomNumber(neighbor);
            if (sharedState.getSourceTerminalSet()[neighbourIndex] < 1) {
                sharedState.getSourceTerminalSet()[neighbourIndex] = size;
                sourceTerminalSize++;
            }
        }

        List<IAtom> targetNeighbours = target.getConnectedAtomsList(target.getAtom(targetAtom));
        for (IAtom neighbor : targetNeighbours) {
            int neighbourIndex = target.getAtomNumber(neighbor);
            if (sharedState.getTargetTerminalSet()[neighbourIndex] < 1) {
                sharedState.getTargetTerminalSet()[neighbourIndex] = size;
                targetTerminalSize++;
            }
        }
    }

    private boolean isEmpty() {
        return candidates.getSourceAtom() == -1 ? true : false;
    }

// Restores the shared state to how it was before adding the last
// candidate pair. Assumes addPair() has been called on the state only once.
    @Override
    public void backTrack() {

        if (isEmpty() || isGoal()) {
            return;
        }
        int addedSourceAtom = candidates.getSourceAtom();

        if (sharedState.getSourceTerminalSet()[addedSourceAtom] == size) {
            sharedState.getSourceTerminalSet()[addedSourceAtom] = 0;
        }

        List<IAtom> sourceNeighbours =
                source.getConnectedAtomsList(source.getAtom(addedSourceAtom));
        for (IAtom neighbor : sourceNeighbours) {
            int neighbourIndex = source.getAtomNumber(neighbor);
            if (sharedState.getSourceTerminalSet()[neighbourIndex] == size) {
                sharedState.getSourceTerminalSet()[neighbourIndex] = 0;
            }
        }

        int addedTargetAtom = candidates.getTargetAtom();

        if (sharedState.getTargetTerminalSet()[addedTargetAtom] == size) {
            sharedState.getTargetTerminalSet()[addedTargetAtom] = 0;
        }

        List<IAtom> targetNeighbours =
                target.getConnectedAtomsList(target.getAtom(addedTargetAtom));
        for (IAtom neighbor : targetNeighbours) {
            int neighbourIndex = target.getAtomNumber(neighbor);
            if (sharedState.getTargetTerminalSet()[neighbourIndex] == size) {
                sharedState.getTargetTerminalSet()[neighbourIndex] = 0;
            }
        }

        sharedState.getSourceMapping()[addedSourceAtom] = -1;
        sharedState.getTargetMapping()[addedTargetAtom] = -1;
        size--;
        candidates = new Match<Integer, Integer>(-1, -1);
    }

    @Override
    public boolean isMatchFeasible(Match<Integer, Integer> candidate) {
        int sourceAtom = candidate.getSourceAtom();
        int targetAtom = candidate.getTargetAtom();

//        int sourceAtomLabel =
//                Integer.parseInt(source.getAtom(sourceAtom).getID());
//        int targetAtomLabel =
//                Integer.parseInt(target.getAtom(targetAtom).getID());
//
//        if (sourceAtomLabel != targetAtomLabel) {
//            return false;
//        }

        if (!matchAtoms(source.getAtom(sourceAtom), target.getAtom(targetAtom))) {
            return false;
        }

        int sourceTerminalNeighborCount = 0;
        int targetTerminalNeighborCount = 0;
        int sourceNewNeighborCount = 0;
        int targetNewNeighborCount = 0;

        List<IAtom> sourceNeighbours =
                source.getConnectedAtomsList(source.getAtom(sourceAtom));
        for (IAtom neighbour : sourceNeighbours) {
            int neighbourIndex = source.getAtomNumber(neighbour);

            IAtom sourceAtomAtom = source.getAtom(sourceAtom);
            IBond sourceBond = source.getBond(sourceAtomAtom, neighbour);

            if (sharedState.getSourceMapping()[neighbourIndex] != -1) {
                int targetNeighbor = sharedState.getSourceMapping()[neighbourIndex];
                IAtom targetNeighbourAtom = target.getAtom(targetNeighbor);
                IAtom targetAtomAtom = target.getAtom(targetAtom);

                if (target.getBond(targetAtomAtom, targetNeighbourAtom) == null) {
                    return false;
                }

                IBond targetBond = target.getBond(targetAtomAtom, targetNeighbourAtom);
                if (!matchBonds(sourceBond, targetBond)) {
                    return false;
                }

            } else {
                if (sharedState.getSourceTerminalSet()[neighbourIndex] > 0) {
                    sourceTerminalNeighborCount++;
                } else {
                    sourceNewNeighborCount++;
                }
            }
        }

        List<IAtom> targetNeighbours =
                target.getConnectedAtomsList(target.getAtom(targetAtom));
        for (IAtom neighbour : targetNeighbours) {
            int neighbourIndex = target.getAtomNumber(neighbour);
            if (sharedState.getTargetMapping()[neighbourIndex] != -1) {
                int sourceNeighbor = sharedState.getTargetMapping()[neighbourIndex];
                IAtom sourceNeighbourAtom = source.getAtom(sourceNeighbor);
                IAtom sourceAtomAtom = source.getAtom(sourceAtom);

                if (source.getBond(sourceAtomAtom, sourceNeighbourAtom) == null) {
                    return false;
                }
            } else {
                if (sharedState.getTargetTerminalSet()[neighbourIndex] > 0) {
                    targetTerminalNeighborCount++;
                } else {
                    targetNewNeighborCount++;
                }
            }
        }
        return (sourceTerminalNeighborCount <= targetTerminalNeighborCount)
                && (sourceNewNeighborCount <= targetNewNeighborCount);
    }

    boolean matchBonds(IBond sourceBond, IBond targetBond) {
        BondMatcher bondMatcher = new VFBondMatcher(source, sourceBond, true);
        return bondMatcher.matches(target, targetBond);
    }

    boolean matchAtoms(IAtom sourceAtom, IAtom targetAtom) {
        AtomMatcher atomMatcher = new VFAtomMatcher(source, sourceAtom, true);
        return atomMatcher.matches(target, targetAtom);
    }
}
