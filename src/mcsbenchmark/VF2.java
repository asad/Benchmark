package mcsbenchmark;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.labelling.AtomContainerPrinter;

public class VF2 {

    // The SharedState class holds four vectors containing the mapping between
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
    }

    class Pair<T, S> {

        T first;
        S second;

        Pair(T a, S b) {
            this.first = a;
            this.second = b;
        }

        @Override
        public String toString() {
            return "(" + first + ", " + second + ")";
        }
    }

    class AtomMapping {

        private IAtomContainer a;
        private IAtomContainer b;
        private Map<IAtom, IAtom> mapping;

        public AtomMapping(IAtomContainer a, IAtomContainer b) {
            this.a = a;
            this.b = b;
            this.mapping = new HashMap<IAtom, IAtom>();
        }

        public void add(IAtom atom1, IAtom atom2) {
            mapping.put(atom1, atom2);
        }

        @Override
        public String toString() {
            String s = "[";
            for (IAtom key : mapping.keySet()) {
                int keyIndex = a.getAtomNumber(key);
                int valueIndex = b.getAtomNumber(mapping.get(key));
                s += keyIndex + ":" + valueIndex + "|";
            }
            return s + "]";
        }

        public boolean isEmpty() {
            return mapping.isEmpty();
        }
    }

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
        // to the
        // state. The candidate should be checked for feasibility and then added
        // using
        // the addPair() method.
        Pair<Integer, Integer> nextCandidate(
                Pair<Integer, Integer> lastCandidate) {
            int lastSourceAtom = lastCandidate.first;
            int lastTargetAtom = lastCandidate.second;

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
        // candidate
        // pair must be feasible to add it to the state.
        void addPair(Pair<Integer, Integer> candidate) {
            size++;
            lastAddition = candidate;

            int sourceAtom = candidate.first;
            int targetAtom = candidate.second;

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
        // candidate
        // pair. Assumes addPair() has been called on the state only once.
        void backTrack() {
//            System.out.println("backtracking " + lastAddition);
            if (lastAddition.first == -1) {
                return;   // XXX hack
            }
            int addedSourceAtom = lastAddition.first;

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

            int addedTargetAtom = lastAddition.second;

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
            int sourceAtom = candidate.first;
            int targetAtom = candidate.second;

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

//            System.out.println("Source : " + Arrays.toString(sharedState.sourceMapping));
//            System.out.println("Target : " + Arrays.toString(sharedState.targetMapping));
//            System.out.println("Source atom " + sourceAtom);
            List<IAtom> sourceNeighbours =
                    source.getConnectedAtomsList(source.getAtom(sourceAtom));
            for (IAtom neighbour : sourceNeighbours) {
                int neighbourIndex = source.getAtomNumber(neighbour);
//                System.out.println("Source neighbour " + neighbourIndex);

                IAtom sourceAtomAtom = source.getAtom(sourceAtom);
                IBond sourceBond = source.getBond(sourceAtomAtom, neighbour);

                if (sharedState.sourceMapping[neighbourIndex] != -1) {
                    int targetNeighbor = sharedState.sourceMapping[neighbourIndex];
//                    System.out.println("targetNeighbour = " + targetNeighbor);
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
//                    System.out.println("Not mapped sourceTerminalSet = " + sharedState.sourceTerminalSet[neighbourIndex]);
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

//            System.out.println(sourceTerminalNeighborCount + " "
//                    + targetTerminalNeighborCount + " "
//                    + sourceNewNeighborCount + " "
//                    + targetNewNeighborCount);
            return (sourceTerminalNeighborCount <= targetTerminalNeighborCount)
                    && (sourceNewNeighborCount <= targetNewNeighborCount);
        }
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

            if (candidate.first == -1) {
                return false;
            }

            lastCandidate = candidate;
//            System.out.println("lastCandidate " + lastCandidate);

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

    void setIDs(IAtomContainer atomContainer) {
        for (int i = 0; i < atomContainer.getAtomCount(); i++) {
            atomContainer.getAtom(i).setID(String.valueOf(i));
        }
    }

    // The isomorphism_vf2() method returns an isomorphism between two molecular
    // graphs using the VF2 algorithm. This can be used for finding both
    // graph-graph isomorphisms and graph-subgraph isomorphisms. In the latter
    // case
    // graph 'a' is the subgraph, implying a.size() < b.size(). In the case that
    // no isomorphism is found an empty mapping is returned.
    AtomMapping isomorphism(IAtomContainer a, IAtomContainer b) {
//        AtomContainerPrinter printer = new AtomContainerPrinter();
        setIDs(a);
//        System.out.println(printer.toString(a));
        setIDs(b);
//        System.out.println(printer.toString(b));
        List<AtomMapping> mappings = new ArrayList<AtomMapping>();
        if (testIsSubgraphHeuristics(a, b)) {
            State state = new State(a, b);
            match(state, mappings);
        }
//        System.out.println("mapping count " + mappings.size());
        return mappings.isEmpty() ? new AtomMapping(a, b) : mappings.get(0);
    }

    /**
     *  Checks some simple heuristics for whether the subgraph query can
     *  realistically be atom subgraph of the supergraph. If, for example, the
     *  number of nitrogen atoms in the query is larger than that of the supergraph
     *  it cannot be part of it.
     *
     * @param  ac1  the supergraph to be checked. 
     * @param  ac2  the subgraph to be tested for. Must not be an IQueryAtomContainer.
     * @return    true if the subgraph ac1 has atom chance to be atom subgraph of ac2
     * @throws org.openscience.cdk.exception.CDKException if the first molecule is an instance
     * of IQueryAtomContainer
     */
    private static boolean testIsSubgraphHeuristics(IAtomContainer ac1, IAtomContainer ac2) {

        int ac1SingleBondCount = 0;
        int ac1DoubleBondCount = 0;
        int ac1TripleBondCount = 0;
        int ac1AromaticBondCount = 0;
        int ac2SingleBondCount = 0;
        int ac2DoubleBondCount = 0;
        int ac2TripleBondCount = 0;
        int ac2AromaticBondCount = 0;

        IBond bond = null;

        for (int i = 0; i < ac1.getBondCount(); i++) {
            bond = ac1.getBond(i);
            if (bond.getFlag(CDKConstants.ISAROMATIC)) {
                ac1AromaticBondCount++;
            } else if (bond.getOrder() == IBond.Order.SINGLE) {
                ac1SingleBondCount++;
            } else if (bond.getOrder() == IBond.Order.DOUBLE) {
                ac1DoubleBondCount++;
            } else if (bond.getOrder() == IBond.Order.TRIPLE) {
                ac1TripleBondCount++;
            }
        }
        for (int i = 0; i < ac2.getBondCount(); i++) {
            bond = ac2.getBond(i);

            if (bond.getFlag(CDKConstants.ISAROMATIC)) {
                ac2AromaticBondCount++;
            } else if (bond.getOrder() == IBond.Order.SINGLE) {
                ac2SingleBondCount++;
            } else if (bond.getOrder() == IBond.Order.DOUBLE) {
                ac2DoubleBondCount++;
            } else if (bond.getOrder() == IBond.Order.TRIPLE) {
                ac2TripleBondCount++;
            }
        }

        if (ac2SingleBondCount < ac1SingleBondCount) {
            return false;
        }
        if (ac2AromaticBondCount < ac1AromaticBondCount) {
            return false;
        }
        if (ac2DoubleBondCount < ac1DoubleBondCount) {
            return false;
        }
        if (ac2TripleBondCount < ac1TripleBondCount) {
            return false;
        }

        IAtom atom = null;
        Map<String, Integer> map = new HashMap<String, Integer>();
        for (int i = 0; i < ac1.getAtomCount(); i++) {
            atom = ac1.getAtom(i);

            if (map.containsKey(atom.getSymbol())) {
                int val = map.get(atom.getSymbol()) + 1;
                map.put(atom.getSymbol(), val);
            } else {
                map.put(atom.getSymbol(), 1);
            }
        }
        for (int i = 0; i < ac2.getAtomCount(); i++) {
            atom = ac2.getAtom(i);
            if (map.containsKey(atom.getSymbol())) {
                int val = map.get(atom.getSymbol()) - 1;
                if (val > 0) {
                    map.put(atom.getSymbol(), val);
                } else {
                    map.remove(atom.getSymbol());
                }
            }
        }
//        System.out.println("Map " + map);
        return map.isEmpty();
    }

    public static void main(String[] args) throws InvalidSmilesException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CCOC(=O)c1[nH]c2ccc(C)cc2(c1(N))");
        IAtomContainer target = sp.parseSmiles("CCOC(=O)c2[nH]c1ccc(C)cc1c2(N=CN(CC)CC)");
////        IAtomContainer target = sp.parseSmiles("CCOC(=O)c1[nH]c2ccc(C)cc2(c1(N))");
//        IAtomContainer query = sp.parseSmiles("CCO");
//        IAtomContainer target = sp.parseSmiles("CCOCN");
        if (query.getAtomCount() <= target.getAtomCount()) {
            VF2 matcher = new VF2();
            AtomMapping mapping = matcher.isomorphism(query, target);
            System.out.println("mapping " + mapping);
        }
    }
}
