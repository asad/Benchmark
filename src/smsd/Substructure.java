/* Copyright (C) 2009-2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
package smsd;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import smsd.algorithm.single.SingleMappingHandler;
import smsd.algorithm.vflib.VF2Sub;
import chemkit.substructure.VF2;

/**
 * This is an ultra fast method to report if query
 * is a substructure for target molecule. If this case is true
 * then it returns only all mapping.
 *
 * This is much faster than {@link
 * org.openscience.cdk.smsd.algorithm.vflib.substructure} class
 * as it only reports first match and backtracks.
 *
 * This class should only be used to report if a query
 * graph is a substructure of the target graph.
 * 
 *  *
 * <p>An example for <b>Substructure search</b>:</p>
 *  <font color="#003366">
 *  <pre>
 * SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
 * IAtomContainer query = sp.parseSmiles("CC");
 * IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
 * Substructure smsd = new Substructure(query, target, true, false, true);
 * Assert.assertTrue(smsd.isSubgraph());
 * Assert.assertEquals(18, smsd.getAllAtomMapping().size());

 * IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
 * smsd = new Substructure(queryContainer, target, true);
 * Assert.assertTrue(smsd.isSubgraph());
 *
 *  </pre>
 *  </font>
 *
 *
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public final class Substructure extends BaseMapping {

    private int vfMappingSize = -1;
    private final ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(Substructure.class);

    /**
     * Constructor for VF Substructure Algorithm 
     * @param query
     * @param target
     * @param shouldMatchBonds Match bond types (i.e. double to double etc)
     * @param matchRings Match ring atoms and ring size
     * @param findAllSubgraph report all subgraphs
     * @throws CDKException  
     */
    public Substructure(
            IAtomContainer query,
            IAtomContainer target,
            boolean shouldMatchBonds,
            boolean matchRings,
            boolean findAllSubgraph) throws CDKException {
        this.mol1 = query;
        this.mol2 = target;
        this.mcsList = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        this.matchBonds = shouldMatchBonds;
        this.matchRings = matchRings;
        if (findAllSubgraph) {
            this.subgraph = findSubgraphs();
        } else {
            this.subgraph = findSubgraph();
        }
    }

    /**
     * Constructor for VF Substructure Algorithm 
     * @param query
     * @param target
     * @param findAllSubgraph report all subgraphs
     * @throws CDKException  
     */
    public Substructure(
            IQueryAtomContainer query,
            IAtomContainer target,
            boolean findAllSubgraph) throws CDKException {
        this.mol1 = query;
        this.mol2 = target;
        this.matchBonds = true;
        this.matchRings = true;
        this.mcsList = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        if (findAllSubgraph) {
            this.subgraph = findSubgraphs();
        } else {
            this.subgraph = findSubgraph();
        }
    }

    private synchronized boolean hasMap(AtomAtomMapping map, List<AtomAtomMapping> mapGlobal) {
        for (AtomAtomMapping test : mapGlobal) {
            if (test.equals(map)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns true if query is a subgraph of target molecule
     * @return
     * @throws CDKException  
     */
    private synchronized boolean findSubgraph() throws CDKException {
        boolean isSubgraph = false;

        if ((mol2 == null) || (mol1 == null)) {
            throw new CDKException("Query or Target molecule is not initialized (NULL)");
        }

        if (mol1.getAtomCount() == 1 || mol2.getAtomCount() == 1) {
            singleMapping(isMatchBonds());
        } else {
            if (mol1.getAtomCount() > mol2.getAtomCount()) {
                return false;
            }
            VF2 mapper = new VF2(isMatchBonds(), isMatchRings());
            List<AtomAtomMapping> mappingsVF2 = new ArrayList<AtomAtomMapping>();
            if (mol1 instanceof IQueryAtomContainer) {
                mapper.set((IQueryAtomContainer) mol1, mol2);
            } else {
                mapper.set(mol1, mol2);
            }
            isSubgraph = mapper.isSubgraph();
            List<AtomAtomMapping> atomMappings = mapper.getAllAtomMapping();
//            System.out.println("Mapping Size " + atomMappings.getCount());
            if (isSubgraph) {
                mappingsVF2.addAll(atomMappings);
            } else {
                return false;
            }
            setVFMappings(mappingsVF2);
        }
        return isSubgraph;
    }

    /**
     * Returns true if query is a subgraph of target molecule
     * @return
     * @throws CDKException  
     */
    private synchronized boolean findSubgraphs() throws CDKException {
        boolean isSubgraph = false;

        if ((mol2 == null) || (mol1 == null)) {
            throw new CDKException("Query or Target molecule is not initialized (NULL)");
        }

        if (mol1.getAtomCount() == 1 || mol2.getAtomCount() == 1) {
            singleMapping(isMatchBonds());
        } else {
            if (mol1.getAtomCount() > mol2.getAtomCount()) {
                return false;
            } else {
                VF2Sub mapper = new VF2Sub(isMatchBonds(), isMatchRings());
                List<AtomAtomMapping> mappingsVF2 = new ArrayList<AtomAtomMapping>();
                if (mol1 instanceof IQueryAtomContainer) {
                    mapper.set((IQueryAtomContainer) mol1, mol2);
                } else {
                    mapper.set(mol1, mol2);
                }
                isSubgraph = mapper.isSubgraph();
                List<AtomAtomMapping> atomMappings = mapper.getAllAtomMapping();
//                    System.out.println("Mapping Size " + atomMapping.getCount());
                if (isSubgraph) {
                    mappingsVF2.addAll(atomMappings);
                } else {
                    return false;
                }
                setVFMappings(mappingsVF2);
            }
        }
        return isSubgraph;
    }

    private synchronized void setVFMappings(List<AtomAtomMapping> mappingsVF2) {
        int counter = 0;
        for (AtomAtomMapping solution : mappingsVF2) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(mol1, mol2);
            if (solution.getCount() > vfMappingSize) {
                this.vfMappingSize = solution.getCount();
                counter = 0;
            }
            for (Map.Entry<IAtom, IAtom> mapping : solution.getMappings().entrySet()) {
                IAtom qAtom = null;
                IAtom tAtom = null;

                qAtom = mapping.getKey();
                tAtom = mapping.getValue();

                if (qAtom != null && tAtom != null) {
                    atomatomMapping.put(qAtom, tAtom);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to NULL");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            }
            if (!atomatomMapping.isEmpty() && !hasMap(atomatomMapping, mcsList)
                    && atomatomMapping.getCount() == vfMappingSize) {
                mcsList.add(counter, atomatomMapping);
                counter++;
            }
        }
    }

    private synchronized void singleMapping(boolean shouldMatchBonds) {
        SingleMappingHandler mcs = null;
        mcs = new SingleMappingHandler();
        mcs.set(mol1, mol2);
        mcs.searchMCS(shouldMatchBonds, isMatchRings());
        mcsList.addAll(mcs.getAllAtomMapping());
    }
}
