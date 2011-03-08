/* Copyright (C) 2009-2010  Syed Asad Rahman <asad@ebi.ac.uk>
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
 */
package mcsbenchmark;

import java.util.ArrayList;
import java.util.Collections;
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
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.rgraph.CDKMCS;
import org.openscience.smsd.algorithm.vflib.interfaces.IMapper;
import org.openscience.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.smsd.algorithm.vflib.map.VFMapper;
import org.openscience.smsd.algorithm.vflib.query.QueryCompiler;
import org.openscience.smsd.global.TimeOut;

/**
 * This is an ultra fast method to report if query
 * is a substructure for target molecule. If this case is true
 * then it returns only all mapping.
 *
 * This is much faster than {@link
 * org.openscience.cdk.smsd.algorithm.vflib.VFlibHandler} class
 * as it only reports first match and backtracks.
 *
 * This class should only be used to report if a query
 * graph is a substructure of the target graph.
 *
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Substructure {

    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static Map<Integer, Integer> firstMCS = null;
    private static List<Map<Integer, Integer>> allMCS = null;
    private IQueryAtomContainer queryMol = null;
    private IAtomContainer mol1 = null;
    private IAtomContainer mol2 = null;
    private List<Map<INode, IAtom>> vfLibSolutions = null;
    private int vfMCSSize = -1;
    private boolean bond_Match_Flag = false;
    private final static ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(Substructure.class);

    /**
     * Constructor for VF Substructure Algorithm 
     */
    public Substructure() {
        allAtomMCS = new ArrayList<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new ArrayList<Map<Integer, Integer>>();

        TimeOut tmo = TimeOut.getInstance();
        tmo.setCDKMCSTimeOut(0.15);
        tmo.setMCSPlusTimeout(0.15);
        tmo.setVFTimeout(0.15);
    }

    private void setFirstMappings() {
        if (!allAtomMCS.isEmpty()) {
            atomsMCS.putAll(allAtomMCS.get(0));
            firstMCS.putAll(allMCS.get(0));
        }
    }

    /** {@inheritDoc}
     *
     * Set the VFLib MCS software
     *
     * @param reactant
     * @param product
     */
    public void set(IAtomContainer reactant, IAtomContainer product) {
        this.mol1 = reactant;
        this.mol2 = product;
    }

    /** {@inheritDoc}
     *
     * @param source
     * @param target
     */
    public void set(IQueryAtomContainer source, IAtomContainer target) {
        queryMol = source;
        mol2 = target;
    }

    private boolean hasMap(Map<Integer, Integer> map, List<Map<Integer, Integer>> mapGlobal) {
        for (Map<Integer, Integer> test : mapGlobal) {
            if (test.equals(map)) {
                return true;
            }
        }
        return false;
    }

    /** {@inheritDoc}
     *
     * @return 
     */
    public List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return Collections.unmodifiableList(allAtomMCS);
    }

    /** {@inheritDoc}
     * @return 
     */
    public List<Map<Integer, Integer>> getAllMapping() {
        return Collections.unmodifiableList(allMCS);
    }

    /** {@inheritDoc}
     * @return 
     */
    public Map<IAtom, IAtom> getFirstAtomMapping() {
        return Collections.unmodifiableMap(atomsMCS);
    }

    /** {@inheritDoc}
     * @return 
     */
    public Map<Integer, Integer> getFirstMapping() {
        return Collections.unmodifiableMap(firstMCS);
    }

    private void setVFMappings(IQuery query) {
        int counter = 0;
        for (Map<INode, IAtom> solution : vfLibSolutions) {
            Map<IAtom, IAtom> atomatomMapping = new HashMap<IAtom, IAtom>();
            Map<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();
            if (solution.size() > vfMCSSize) {
                this.vfMCSSize = solution.size();
                counter = 0;
            }
            for (Map.Entry<INode, IAtom> mapping : solution.entrySet()) {
                IAtom qAtom = null;
                IAtom tAtom = null;

                qAtom = query.getAtom(mapping.getKey());
                tAtom = mapping.getValue();

                Integer qIndex = Integer.valueOf(getReactantMol().getAtomNumber(qAtom));
                Integer tIndex = Integer.valueOf(getProductMol().getAtomNumber(tAtom));
                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to NULL");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            }
            if (!atomatomMapping.isEmpty() && !hasMap(indexindexMapping, allMCS)
                    && indexindexMapping.size() == vfMCSSize) {
                allAtomMCS.add(counter, atomatomMapping);
                allMCS.add(counter, indexindexMapping);
                counter++;
            }
        }
    }

    public boolean isSubgraph(boolean shouldMatchBonds) throws CDKException {

        setBondMatchFlag(shouldMatchBonds);
        if (getReactantMol().getAtomCount() > getProductMol().getAtomCount()) {
            return false;
        } else if ((isBondMatchFlag() && testIsSubgraphHeuristics(getReactantMol(), getProductMol()))
                || !isBondMatchFlag()) {
            //boolean flag = CDKMCS.isSubgraph(getProductMol(), getReactantMol(), shouldMatchBonds);
//            if (!CDKMCS.isTimeOut()) {
//                return flag;
//            } else {
                IQuery query = null;
                IMapper mapper = null;
                vfLibSolutions = new ArrayList<Map<INode, IAtom>>();
                if (queryMol != null) {
                    query = new QueryCompiler(queryMol).compile();
                    mapper = new VFMapper(query);
                    if (mapper.hasMap(getProductMol())) {
                        return true;
                    }
                } else {
                    query = new QueryCompiler(mol1, isBondMatchFlag()).compile();
                    mapper = new VFMapper(query);
                    if (mapper.hasMap(getProductMol())) {
                        return true;
                    }
                }
                return false;
//            }
        } else {
            return false;
        }
    }

    /**
     * 
     * @param shouldMatchBonds
     * @return
     */
    public boolean findSubgraphs(boolean shouldMatchBonds) {

        setBondMatchFlag(shouldMatchBonds);
        if (getReactantMol().getAtomCount() > getProductMol().getAtomCount()) {
            return false;
        } else if ((isBondMatchFlag() && testIsSubgraphHeuristics(getReactantMol(), getProductMol()))
                || !isBondMatchFlag()) {

            IQuery query = null;
            IMapper mapper = null;
            vfLibSolutions = new ArrayList<Map<INode, IAtom>>();
            if (queryMol != null) {
                query = new QueryCompiler(queryMol).compile();
                mapper = new VFMapper(query);
                if (mapper.hasMap(getProductMol())) {
                    List<Map<INode, IAtom>> maps = mapper.getMaps(getProductMol());
                    if (maps != null) {
                        vfLibSolutions.addAll(maps);
                    }
                } else {
                    return false;
                }

            } else {
                query = new QueryCompiler(mol1, isBondMatchFlag()).compile();
                mapper = new VFMapper(query);
                if (mapper.hasMap(getProductMol())) {
                    List<Map<INode, IAtom>> maps = mapper.getMaps(getProductMol());
                    if (maps != null) {
                        vfLibSolutions.addAll(maps);
                    }
                } else {
                    return false;
                }
            }
            setVFMappings(query);
        }
        if (!allAtomMCS.isEmpty()) {
            setFirstMappings();
        }
        return (!allMCS.isEmpty() && allMCS.iterator().next().size() == getReactantMol().getAtomCount()) ? true : false;
    }

    /**
     * 
     * @param shouldMatchBonds
     * @return
     */
    public boolean findSubgraph(boolean shouldMatchBonds) {

        setBondMatchFlag(shouldMatchBonds);
        if (getReactantMol().getAtomCount() > getProductMol().getAtomCount()) {
            return false;
        } else if ((isBondMatchFlag() && testIsSubgraphHeuristics(getReactantMol(), getProductMol()))
                || !isBondMatchFlag()) {

            IQuery query = null;
            IMapper mapper = null;
            vfLibSolutions = new ArrayList<Map<INode, IAtom>>();
            if (queryMol != null) {
                query = new QueryCompiler(queryMol).compile();
                mapper = new VFMapper(query);
                if (mapper.hasMap(getProductMol())) {
                    Map<INode, IAtom> map = mapper.getFirstMap(getProductMol());
                    if (map != null) {
                        vfLibSolutions.add(map);
                    }
                } else {
                    return false;
                }

            } else {
                query = new QueryCompiler(mol1, isBondMatchFlag()).compile();
                mapper = new VFMapper(query);
                if (mapper.hasMap(getProductMol())) {
                    Map<INode, IAtom> map = mapper.getFirstMap(getProductMol());
                    if (map != null) {
                        vfLibSolutions.add(map);
                    }
                } else {
                    return false;
                }
            }
            setVFMappings(query);
        }
        if (!allAtomMCS.isEmpty()) {
            setFirstMappings();
        }
        return (!allMCS.isEmpty() && allMCS.iterator().next().size() == getReactantMol().getAtomCount()) ? true : false;
    }

    /**
     * @return the shouldMatchBonds
     */
    public boolean isBondMatchFlag() {
        return bond_Match_Flag;
    }

    /**
     * @param shouldMatchBonds the shouldMatchBonds to set
     */
    public void setBondMatchFlag(boolean shouldMatchBonds) {
        this.bond_Match_Flag = shouldMatchBonds;
    }

    private IAtomContainer getReactantMol() {
        return queryMol == null ? mol1 : queryMol;
    }

    private IAtomContainer getProductMol() {
        return mol2;
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
}
