package vf2;

import java.util.List;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

// The State class represents a single state in the isomorphism detection
// algorithm. Every state uses and modifies the same SharedState object.
class State {

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
    Pair<Integer, Integer> lastAddition;
    SharedState sharedState;
    boolean ownSharedState;

    State(IAtomContainer source, IAtomContainer target) {
        this.size = 0;
        this.sourceTerminalSize = 0;
        this.targetTerminalSize = 0;
        this.source = source;
        this.target = target;
        this.lastAddition = new Pair<Integer, Integer>(-1, -1);
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
        this.lastAddition = new Pair<Integer, Integer>(-1, -1);
        this.sharedState = state.sharedState;
        this.ownSharedState = false;
    }

    // Returns true if the state contains an isomorphism.
    boolean succeeded() {
        return size == source.getAtomCount();
    }

    // Returns the current isomorphism for the state in an AtomMapping
    // object.
    AtomMapping getMapping() {
        AtomMapping mapping = new AtomMapping(source, target);

        for (int i = 0; i < size; i++) {
            mapping.add(source.getAtom(i),
                    target.getAtom(sharedState.sourceMapping[i]));
        }
        return mapping;
    }

    // Returns the next candidate pair (sourceAtom, targetAtom) to be added
    // to the state. The candidate should be checked for feasibility and then added
    // using the addPair() method.
    Pair<Integer, Integer> nextCandidate(
            Pair<Integer, Integer> lastCandidate) {
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
                    && (sharedState.sourceMapping[lastSourceAtom] != -1 || sharedState.sourceTerminalSet[lastSourceAtom] == 0)) {
                lastSourceAtom++;
                lastTargetAtom = 0;
            }
        } else {
            while (lastSourceAtom < sourceSize
                    && sharedState.sourceMapping[lastSourceAtom] != -1) {
                lastSourceAtom++;
                lastTargetAtom = 0;
            }
        }

        if (sourceTerminalSize > size && targetTerminalSize > size) {
            while (lastTargetAtom < targetSize
                    && (sharedState.targetMapping[lastTargetAtom] != -1 || sharedState.targetTerminalSet[lastTargetAtom] == 0)) {
                lastTargetAtom++;
            }
        } else {
            while (lastTargetAtom < targetSize
                    && sharedState.targetMapping[lastTargetAtom] != -1) {
                lastTargetAtom++;
            }
        }

        if (lastSourceAtom < sourceSize && lastTargetAtom < targetSize) {
            return new Pair<Integer, Integer>(lastSourceAtom, lastTargetAtom);
        }

        return new Pair<Integer, Integer>(-1, -1);
    }

    // Adds the candidate pair (sourceAtom, targetAtom) to the state. The
    // candidate pair must be feasible to add it to the state.
    void addPair(Pair<Integer, Integer> candidate) {
        size++;
        lastAddition = candidate;

        int sourceAtom = candidate.getSourceAtom();
        int targetAtom = candidate.getTargetAtom();

        if (sharedState.sourceTerminalSet[sourceAtom] < 1) {
            sharedState.sourceTerminalSet[sourceAtom] = size;
//                sourceTerminalSize++;
        }

        if (sharedState.targetTerminalSet[targetAtom] < 1) {
            sharedState.targetTerminalSet[targetAtom] = size;
//                targetTerminalSize++;
        }

        sharedState.sourceMapping[sourceAtom] = targetAtom;
        sharedState.targetMapping[targetAtom] = sourceAtom;

        List<IAtom> sourceNeighbours =
                source.getConnectedAtomsList(source.getAtom(sourceAtom));
        for (IAtom neighbor : sourceNeighbours) {
            int neighbourIndex = source.getAtomNumber(neighbor);
            if (sharedState.sourceTerminalSet[neighbourIndex] < 1) {
                sharedState.sourceTerminalSet[neighbourIndex] = size;
                sourceTerminalSize++;
            }
        }

        List<IAtom> targetNeighbours = target.getConnectedAtomsList(target.getAtom(targetAtom));
        for (IAtom neighbor : targetNeighbours) {
            int neighbourIndex = target.getAtomNumber(neighbor);
            if (sharedState.targetTerminalSet[neighbourIndex] < 1) {
                sharedState.targetTerminalSet[neighbourIndex] = size;
                targetTerminalSize++;
            }
        }
    }

    // Restores the shared state to how it was before adding the last
    // candidate pair. Assumes addPair() has been called on the state only once.
    void backTrack() {
        if (lastAddition.getSourceAtom() == -1) {
            return;   // XXX hack
        }
        int addedSourceAtom = lastAddition.getSourceAtom();

        if (sharedState.sourceTerminalSet[addedSourceAtom] == size) {
            sharedState.sourceTerminalSet[addedSourceAtom] = 0;
        }

        List<IAtom> sourceNeighbours =
                source.getConnectedAtomsList(source.getAtom(addedSourceAtom));
        for (IAtom neighbor : sourceNeighbours) {
            int neighbourIndex = source.getAtomNumber(neighbor);
            if (sharedState.sourceTerminalSet[neighbourIndex] == size) {
                sharedState.sourceTerminalSet[neighbourIndex] = 0;
            }
        }

        int addedTargetAtom = lastAddition.getTargetAtom();

        if (sharedState.targetTerminalSet[addedTargetAtom] == size) {
            sharedState.targetTerminalSet[addedTargetAtom] = 0;
        }

        List<IAtom> targetNeighbours =
                target.getConnectedAtomsList(target.getAtom(addedTargetAtom));
        for (IAtom neighbor : targetNeighbours) {
            int neighbourIndex = target.getAtomNumber(neighbor);
            if (sharedState.targetTerminalSet[neighbourIndex] == size) {
                sharedState.targetTerminalSet[neighbourIndex] = 0;
            }
        }

        sharedState.sourceMapping[addedSourceAtom] = -1;
        sharedState.targetMapping[addedTargetAtom] = -1;
        size--;
        lastAddition = new Pair<Integer, Integer>(-1, -1);
    }

    boolean isFeasible(Pair<Integer, Integer> candidate) {
        int sourceAtom = candidate.getSourceAtom();
        int targetAtom = candidate.getTargetAtom();

        int sourceAtomLabel =
                Integer.parseInt(source.getAtom(sourceAtom).getID());
        int targetAtomLabel =
                Integer.parseInt(target.getAtom(targetAtom).getID());

        if (sourceAtomLabel != targetAtomLabel) {
            return false;
        }

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

            if (sharedState.sourceMapping[neighbourIndex] != -1) {
                int targetNeighbor = sharedState.sourceMapping[neighbourIndex];
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
                if (sharedState.sourceTerminalSet[neighbourIndex] > 0) {
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
            if (sharedState.targetMapping[neighbourIndex] != -1) {
                int sourceNeighbor = sharedState.targetMapping[neighbourIndex];
                IAtom sourceNeighbourAtom = source.getAtom(sourceNeighbor);
                IAtom sourceAtomAtom = source.getAtom(targetAtom);

                if (source.getBond(sourceAtomAtom, sourceNeighbourAtom) == null) {
                    return false;
                }
            } else {
                if (sharedState.targetTerminalSet[neighbourIndex] > 0) {
                    targetTerminalNeighborCount++;
                } else {
                    targetNewNeighborCount++;
                }
            }
        }
        return (sourceTerminalNeighborCount <= targetTerminalNeighborCount)
                && (sourceNewNeighborCount <= targetNewNeighborCount);
    }

    boolean match(State state, List<AtomMapping> mappings) {
//        System.out.println("Matched " + state.size + " out of " + state.source.getAtomCount());
        if (state.succeeded()) {
            mappings.add(state.getMapping());
            return true;
        }

        Pair<Integer, Integer> lastCandidate = new Pair<Integer, Integer>(-1, -1);

        boolean found = false;
        while (!found) {
            Pair<Integer, Integer> candidate = state.nextCandidate(lastCandidate);

            if (candidate.getSourceAtom() == -1) {
                return false;
            }

            lastCandidate = candidate;

            if (state.isFeasible(candidate)) {
                State nextState = state;
                nextState.addPair(candidate);
                found = match(nextState, mappings);
                if (found) {
                    return true;
                }
                nextState.backTrack();
            }
        }

        return found;
    }

    boolean matchBonds(IBond sourceBond, IBond targetBond) {
        if ((sourceBond.getFlag(CDKConstants.ISAROMATIC) == targetBond.getFlag(CDKConstants.ISAROMATIC))
                && (sourceBond.getOrder() == targetBond.getOrder())) {
            return true;
        } else if (sourceBond.getFlag(CDKConstants.ISAROMATIC) && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
            return true;
        }

//        System.out.println("Bond order mismatch "
//                + sourceBond.getOrder() + " " + targetBond.getOrder());
        return false;
    }

    boolean matchAtoms(IAtom sourceAtom, IAtom targetAtom) {
        return sourceAtom.getSymbol().equals(targetAtom.getSymbol()) ? true : false;
    }
}
