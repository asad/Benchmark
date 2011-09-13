/*
 *
 *
 * Copyright (C) 2009-2011  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 * 
 ** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
 **
 ** This file is part of chemkit. For more information see
 ** <http://www.chemkit.org>.
 **
 ** chemkit is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU Lesser General Public License as published by
 ** the Free Software Foundation, either version 3 of the License, or
 ** (at your option) any later version.
 **
 ** chemkit is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ** GNU Lesser General Public License for more details.
 **
 ** You should have received a copy of the GNU Lesser General Public License
 ** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
 **
 ******************************************************************************/
package chemkit.substructure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import smsd.AtomAtomMapping;
import smsd.interfaces.AbstractSubGraph;
import smsd.interfaces.IMCSBase;

/**
 * This class finds mapping states between query and target
 * molecules.
 * 
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class VF2 extends AbstractSubGraph implements IMCSBase {

    private List<AtomAtomMapping> allAtomMCS = null;
    private List<Map<Integer, Integer>> allMCS = null;
    private IAtomContainer source = null;
    private IAtomContainer target = null;
    private final boolean shouldMatchRings;
    private final ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(VF2.class);
    private final boolean shouldMatchBonds;

    /**
     * Constructor for an extended VF Algorithm for the MCS search
     * @param shouldMatchBonds
     * @param shouldMatchRings  
     */
    public VF2(boolean shouldMatchBonds, boolean shouldMatchRings) {
        this.shouldMatchRings = shouldMatchRings;
        this.shouldMatchBonds = shouldMatchBonds;
        allAtomMCS = new ArrayList<AtomAtomMapping>();
        allMCS = new ArrayList<Map<Integer, Integer>>();
    }

    /** 
     * The isomorphism method returns an isomorphism between two molecular
     *  graphs using the VF2Automorphism algorithm. This can be used for finding both
     *  graph-graph isomorphisms and graph-subgraph isomorphisms. In the latter
     *  case graph 'a' is the subgraph, implying a.size() < b.size(). In the case that
     *  no isomorphism is found an empty mapping is returned.
     * 
     * 
     * @param shouldMatchBonds 
     * @param shouldMatchRings 
     * @return
     */
    private synchronized void isomorphism() {

        if (shouldMatchRings) {
            try {
                initializeMolecule(source);
                initializeMolecule(target);
            } catch (CDKException ex) {
                Logger.error(Level.SEVERE, null, ex);
            }
        }

        if (!isDead(source, target) && testIsSubgraphHeuristics(source, target, shouldMatchBonds)) {
//            AtomContainerPrinter printer = new AtomContainerPrinter();
//            System.out.println(printer.toString(a));
//            System.out.println(printer.toString(b));
            State state = new State(source, target, shouldMatchBonds, shouldMatchRings);
            if (!state.isDead()) {
                state.matchFirst(state, allAtomMCS);
            }
        }
    }

    /** 
     * The isomorphism method returns an isomorphism between two molecular
     *  graphs using the VF2Automorphism algorithm. This can be used for finding both
     *  graph-graph isomorphisms and graph-subgraph isomorphisms. In the latter
     *  case graph 'a' is the subgraph, implying a.size() < b.size(). In the case that
     *  no isomorphism is found an empty mapping is returned.
     * 
     * 
     */
    private synchronized void isomorphisms() {

        if (shouldMatchRings) {
            try {
                initializeMolecule(source);
                initializeMolecule(target);
            } catch (CDKException ex) {
                Logger.error(Level.SEVERE, null, ex);
            }
        }

        if (!isDead(source, target) && testIsSubgraphHeuristics(source, target, shouldMatchBonds)) {
            State state = new State(source, target, shouldMatchBonds, shouldMatchRings);
            if (!state.isDead()) {
                state.matchAll(state, allAtomMCS);
            }
        }
    }

    // Returns true substructure is bigger than the target
    private synchronized boolean isDead(IAtomContainer a, IAtomContainer b) {
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
    private synchronized boolean testIsSubgraphHeuristics(IAtomContainer ac1, IAtomContainer ac2, boolean shouldMatchBonds) {

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
                if (bond instanceof IQueryBond) {
                    continue;
                }
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
            for (int indexI = 0; indexI < ac2.getBondCount(); indexI++) {
                bond = ac2.getBond(indexI);
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
            if (atom instanceof IQueryAtom) {
                continue;
            }
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
        return map.isEmpty();
    }

    @Override
    public boolean isSubgraph() {
        isomorphism();
        return allAtomMCS.isEmpty() ? false : true;
    }

    public boolean isSubgraphs() {
        isomorphisms();
        return allAtomMCS.isEmpty() ? false : true;
    }

    @Override
    public void set(IAtomContainer source, IAtomContainer target) throws CDKException {
        this.source = source;
        this.target = target;
    }

    @Override
    public void set(IQueryAtomContainer source, IAtomContainer target) throws CDKException {
        this.source = source;
        this.target = target;
    }

    @Override
    public List<AtomAtomMapping> getAllAtomMapping() {
        return allAtomMCS;
    }

    @Override
    public List<Map<Integer, Integer>> getAllMapping() {
        return allMCS;
    }

    @Override
    public AtomAtomMapping getFirstAtomMapping() {
        if (allAtomMCS.iterator().hasNext()) {
            return allAtomMCS.iterator().next();
        }
        return new AtomAtomMapping(source, target);
    }

    @Override
    public Map<Integer, Integer> getFirstMapping() {
        if (allMCS.iterator().hasNext()) {
            return allMCS.iterator().next();
        }
        return new TreeMap<Integer, Integer>();
    }
}
