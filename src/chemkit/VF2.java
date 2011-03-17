/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package chemkit;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Asad
 */
public class VF2 {

    /** The isomorphism method returns an isomorphism between two molecular
     *  graphs using the VF2Automorphism algorithm. This can be used for finding both
     *  graph-graph isomorphisms and graph-subgraph isomorphisms. In the latter
     *  case graph 'a' is the subgraph, implying a.size() < b.size(). In the case that
     *  no isomorphism is found an empty mapping is returned.
    
     * 
     * @param a query molecule
     * @param b target molecule
     * @param shouldMatchBonds 
     * @return
     */
    public AtomMapping isomorphism(IAtomContainer a, IAtomContainer b, boolean shouldMatchBonds) {

        List<AtomMapping> mappings = new ArrayList<AtomMapping>();
        if (!isDead(a, b) && testIsSubgraphHeuristics(a, b, shouldMatchBonds)) {
//            AtomContainerPrinter printer = new AtomContainerPrinter();
//            System.out.println(printer.toString(a));
//            System.out.println(printer.toString(b));
            State state = new State(a, b, shouldMatchBonds);
            if (!state.isDead()) {
                state.matchFirst(state, mappings);
            }
        }
        return mappings.isEmpty() ? new AtomMapping(a, b) : mappings.get(0);
    }

    /**
     * 
     * @param a query molecule
     * @param b target molecule
     * @param shouldMatchBonds 
     * @return
     */
    public List<AtomMapping> isomorphisms(IAtomContainer a, IAtomContainer b, boolean shouldMatchBonds) {

        List<AtomMapping> mappings = new ArrayList<AtomMapping>();
        if (!isDead(a, b) && testIsSubgraphHeuristics(a, b, shouldMatchBonds)) {
////            AtomContainerPrinter printer = new AtomContainerPrinter();
////            System.out.println(printer.toString(a));
//            System.out.println(printer.toString(b));
            State state = new State(a, b, shouldMatchBonds);
            if (!state.isDead()) {
                state.matchAll(state, mappings);
            }
        }
        return mappings;
    }

    // Returns true substructure is bigger than teh target
    public boolean isDead(IAtomContainer a, IAtomContainer b) {
        return a.getAtomCount() > b.getAtomCount();
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
    private static boolean testIsSubgraphHeuristics(IAtomContainer ac1, IAtomContainer ac2, boolean shouldMatchBonds) {

        int ac1SingleBondCount = 0;
        int ac1DoubleBondCount = 0;
        int ac1TripleBondCount = 0;
        int ac1AromaticBondCount = 0;
        int ac2SingleBondCount = 0;
        int ac2DoubleBondCount = 0;
        int ac2TripleBondCount = 0;
        int ac2AromaticBondCount = 0;

        IBond bond = null;

        if (shouldMatchBonds) {
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
}