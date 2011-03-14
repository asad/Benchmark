/*
 *
 *
 * Copyright (C) 2009-2010  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * MX Cheminformatics Tools for Java
 *
 * Copyright (c) 2007-2009 Metamolecular, LLC
 *
 * http://metamolecular.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining targetAtom copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */
package isomorphism.vf2.bond;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;

/**
 * This class finds mapping states between query and target
 * molecules.
 * 
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class VFBondState implements IBondState {

    private List<VFBondMatcher> candidates;
    private IAtomContainer query;
    private IAtomContainer target;
    private List<IBond> queryPath;
    private List<IBond> targetPath;
    private BondMapping map;
    private Map<IBond, List<IBond>> neighbourQueryMap;
    private Map<IBond, List<IBond>> neighbourTargetMap;
    /*Needs AdjacencyMatrix as CDK doesnot store the matches*/
    private boolean[][] atomAdjacencyMatrix;
    private boolean[][] bondAdjacencyMatrix;

    /**
     * initialize the VFState with query and target
     * @param query
     * @param target
     */
    public VFBondState(IAtomContainer query, IAtomContainer target) {
        this.query = query;
        this.target = target;
        this.map = new BondMapping(target, query);
        this.candidates = new ArrayList<VFBondMatcher>();


        this.neighbourQueryMap = new HashMap<IBond, List<IBond>>();
        this.neighbourTargetMap = new HashMap<IBond, List<IBond>>();
//        initialize(query, target);

        if (testIsSubgraphHeuristics(query, target)) {
            this.queryPath = new ArrayList<IBond>(query.getAtomCount());
            this.targetPath = new ArrayList<IBond>(target.getAtomCount());
            this.atomAdjacencyMatrix = new boolean[query.getAtomCount()][target.getAtomCount()];
            this.bondAdjacencyMatrix = new boolean[query.getBondCount()][target.getBondCount()];

            for (IBond bondQ : query.bonds()) {
                for (IBond bondT : target.bonds()) {
                    bondAdjacencyMatrix[Integer.parseInt(bondQ.getID())][Integer.parseInt(bondT.getID())] = false;
                }
            }

            loadRootCandidates();
        }
    }

    private VFBondState(VFBondState state, VFBondMatcher match) {
        this.candidates = new ArrayList<VFBondMatcher>();
        this.queryPath = new ArrayList<IBond>(state.queryPath);
        this.targetPath = new ArrayList<IBond>(state.targetPath);

        this.map = state.map;
        this.query = state.query;
        this.target = state.target;

        this.neighbourQueryMap = state.neighbourQueryMap;
        this.neighbourTargetMap = state.neighbourTargetMap;
        this.atomAdjacencyMatrix = state.atomAdjacencyMatrix;
        this.bondAdjacencyMatrix = state.bondAdjacencyMatrix;

        map.add(match.getQueryBond(), match.getTargetBond());
        queryPath.add(match.getQueryBond());
        targetPath.add(match.getTargetBond());
        loadCandidates(match);
    }

    /** {@inheritDoc}
     */
    @Override
    public void backTrack() {
        if (queryPath.isEmpty() || isGoal()) {
            map.clear();
            return;
        }
        if (isHeadMapped()) {
            return;
        }
        map.clear();
        for (int i = 0; i < queryPath.size() - 1; i++) {
            map.add(queryPath.get(i), targetPath.get(i));
        }
    }

    /** {@inheritDoc}
     */
    @Override
    public Map<IBond, IBond> getMap() {
        return new HashMap<IBond, IBond>(map.getMapping());
    }

    /** {@inheritDoc}
     */
    @Override
    public boolean hasNextCandidate() {
        return !candidates.isEmpty();
    }

    /** {@inheritDoc}
     */
    @Override
    public boolean isDead() {
        return query.getBondCount() > target.getBondCount();
    }

    /** {@inheritDoc}
     */
    @Override
    public boolean isGoal() {
        return map.size() == query.getBondCount();
    }

    /** {@inheritDoc}
     */
    @Override
    public boolean isMatchFeasible(VFBondMatcher match) {
        if (map.containsQueryBond(match.getQueryBond())
                || map.containsTargetBond(match.getTargetBond())) {
            return false;
        }

        if (!isNeighbourFeasible(match)) {
            return false;
        }

        if (!checkBondMatrix(match.getQueryBond(), match.getTargetBond())) {
            return false;
        }

        if (!bondMatcher(match)) {
            return false;
        }

        return true;
    }

    /** {@inheritDoc}
     */
    @Override
    public VFBondMatcher nextCandidate() {
        return candidates.remove(candidates.size() - 1);
    }

    /** {@inheritDoc}
     */
    @Override
    public IBondState nextState(VFBondMatcher match) {
        return new VFBondState(this, match);
    }

    private boolean loadRootCandidates() {
        for (IBond qBond : neighbourQueryMap.keySet()) {
            boolean flag = false;
            for (IBond tBond : neighbourTargetMap.keySet()) {
                VFBondMatcher match = new VFBondMatcher(qBond, tBond);
                if (bondMatchers(match)) {
                    candidates.add(match);
                    bondAdjacencyMatrix[Integer.parseInt(qBond.getID())][Integer.parseInt(tBond.getID())] = true;
                    flag = true;
                }
            }
            if (!flag) {
                candidates.clear();
                return false;
            }
        }

//        System.out.println("Compatibility graph " + candidates.size());
        return true;
    }

//@TODO Asad Check the Neighbour count
    private void loadCandidates(VFBondMatcher lastMatch) {
        List<IBond> queryNeighbors = neighbourQueryMap.get(lastMatch.getQueryBond());
        List<IBond> targetNeighbors = neighbourTargetMap.get(lastMatch.getTargetBond());

        for (IBond queryBond : queryNeighbors) {
            for (IBond targetBond : targetNeighbors) {
                if (checkBondMatrix(queryBond, targetBond)) {
                    VFBondMatcher match = new VFBondMatcher(queryBond, targetBond);
                    if (candidateFeasible(match)) {
//                    System.out.println("map " + map.size());
                        candidates.add(match);
                    }
                }
            }
        }
//        System.out.println("candidates " + candidates.size());
    }

    private boolean candidateFeasible(VFBondMatcher candidate) {
        for (IBond queryAtom : map.queryBonds()) {
            if (queryAtom.equals(candidate.getQueryBond())
                    || map.getMappedTargetBond(queryAtom).equals(candidate.getTargetBond())) {
                return false;
            }
        }
        return true;
    }
    //This function is updated by Asad to include more matches

    private boolean bondMatchers(VFBondMatcher match) {
        if (isNeighbourFeasible(match)) {
            return matchBond(match.getQueryBond(), match.getTargetBond());
        }
        return false;
    }

    private boolean isNeighbourFeasible(VFBondMatcher match) {
        if (neighbourQueryMap.get(match.getQueryBond()).size() > neighbourTargetMap.get(match.getTargetBond()).size()) {
            return false;
        }
        return true;
    }

    private boolean bondMatcher(VFBondMatcher match) {
        if (queryPath.isEmpty()) {
            return true;
        }

        if (!matchBondsToHead(match)) {
            return false;
        }

        for (int i = 0; i < queryPath.size() - 1; i++) {
            IBond queryBond = queryPath.get(i);
            if (!hasCommonAtom(queryBond, match.getQueryBond())) {
                continue;
            }
            IBond targetBond = targetPath.get(i);
            if (!hasCommonAtom(targetBond, match.getTargetBond())) {
                return false;
            }
            if (!checkBondMatrix(queryBond, targetBond)) {
                return false;
            }
        }
        return true;
    }

    private boolean isHeadMapped() {
        IBond head = queryPath.get(queryPath.size() - 1);
        List<IBond> queryHeadNeighbors = neighbourQueryMap.get(head);
        for (IBond neighbor : queryHeadNeighbors) {
            if (!map.containsQueryBond(neighbor)) {
                return false;
            }
        }
        return true;
    }

    private boolean matchBondsToHead(VFBondMatcher match) {
        IBond queryHead = getQueryPathHead();
        IBond targetHead = getTargetPathHead();

        if (hasCommonAtom(queryHead, match.getQueryBond()) && hasCommonAtom(targetHead, match.getTargetBond())) {
//            return checkBondMatrix(match.getQueryBond(), match.getTargetBond());
            if (checkBondMatrix(queryHead, targetHead) && checkBondMatrix(match.getQueryBond(), match.getTargetBond())) {
                return true;
            }
        }
        return false;
    }

    /**
     * Determines if two bonds have at least one atom in common.
     *
     * @param  a  first bond
     * @param  b  second bond
     * @return    the symbol of the common atom or "" if
     *            the 2 bonds have no common atom
     */
    private static boolean hasCommonAtom(IBond a, IBond b) {
        return a.contains(b.getAtom(0)) || a.contains(b.getAtom(1));
    }

    private IBond getQueryPathHead() {
        return queryPath.get(queryPath.size() - 1);
    }

    private IBond getTargetPathHead() {
        return targetPath.get(targetPath.size() - 1);
    }

    private boolean matchBond(IBond queryBond, IBond targetBond) {
        if (queryBond instanceof IQueryBond) {
            return ((IQueryBond) queryBond).matches(targetBond) && matchAtoms(queryBond, targetBond) ? true : false;
        } else if ((queryBond.getFlag(CDKConstants.ISAROMATIC) == targetBond.getFlag(CDKConstants.ISAROMATIC))
                && (queryBond.getOrder() == targetBond.getOrder())) {
            return matchAtoms(queryBond, targetBond);
        } else if (queryBond.getFlag(CDKConstants.ISAROMATIC) && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
            return matchAtoms(queryBond, targetBond);
        }
        return false;
    }

    private boolean matchAtoms(IBond sourceBond, IBond targetBond) {
        if (sourceBond instanceof IQueryBond) {
            return ((IQueryBond) sourceBond).matches(targetBond) ? true : false;
        } else {
            return matchSymbol(sourceBond, targetBond);
        }
    }

    private boolean matchSymbol(IBond queryBond, IBond targetBond) {
        if (( // a1 = a2 && b1 = b2
                queryBond.getAtom(0).getSymbol().equals(targetBond.getAtom(0).getSymbol())
                && queryBond.getAtom(1).getSymbol().equals(targetBond.getAtom(1).getSymbol()))
                || ( // a1 = b2 && b1 = a2
                queryBond.getAtom(0).getSymbol().equals(targetBond.getAtom(1).getSymbol())
                && queryBond.getAtom(1).getSymbol().equals(targetBond.getAtom(0).getSymbol()))) {
            return true;
        }
        return false;
    }

    /**
     *  Checks some simple heuristics for whether the subgraph query can
     *  realistically be atom subgraph of the supergraph. If, for example, the
     *  number of nitrogen atoms in the query is larger than that of the supergraph
     *  it cannot be part of it.
     *
     * @param  query  the subgraph to be checked. 
     * @param  target  the super-graph to be tested for. Must not be an IQueryAtomContainer.
     * @return    true if the subgraph query has atom chance to be atom subgraph of target
     * @throws org.openscience.cdk.exception.CDKException if the first molecule is an instance
     * of IQueryAtomContainer
     */
    private boolean testIsSubgraphHeuristics(IAtomContainer query, IAtomContainer target) {

        int ac1SingleBondCount = 0;
        int ac1DoubleBondCount = 0;
        int ac1TripleBondCount = 0;
        int ac1AromaticBondCount = 0;
        int ac2SingleBondCount = 0;
        int ac2DoubleBondCount = 0;
        int ac2TripleBondCount = 0;
        int ac2AromaticBondCount = 0;

        IBond bond = null;

        for (int i = 0; i < query.getBondCount(); i++) {
            bond = query.getBond(i);
            bond.setID(Integer.toString(i));
            IAtom atom0 = bond.getAtom(0);
            IAtom atom1 = bond.getAtom(1);
            List<IBond> a = new ArrayList<IBond>(query.getConnectedBondsList(atom0));
            List<IBond> b = new ArrayList<IBond>(query.getConnectedBondsList(atom1));
            for (IBond bondNeigh : b) {
                if (!a.contains(bondNeigh)) {
                    a.add(bondNeigh);
                }
            }
            neighbourQueryMap.put(bond, a);
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
        for (int i = 0; i < target.getBondCount(); i++) {
            bond = target.getBond(i);
            bond.setID(Integer.toString(i));
            IAtom atom0 = bond.getAtom(0);
            IAtom atom1 = bond.getAtom(1);
            List<IBond> a = new ArrayList<IBond>(target.getConnectedBondsList(atom0));
            List<IBond> b = new ArrayList<IBond>(target.getConnectedBondsList(atom1));
            for (IBond bondNeigh : b) {
                if (!a.contains(bondNeigh)) {
                    a.add(bondNeigh);
                }
            }
            neighbourTargetMap.put(bond, a);

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
        Map<String, Integer> symbolMap = new HashMap<String, Integer>();
        for (int i = 0; i < query.getAtomCount(); i++) {
            atom = query.getAtom(i);
            atom.setID(Integer.toString(i));
            if (symbolMap.containsKey(atom.getSymbol())) {
                int val = symbolMap.get(atom.getSymbol()) + 1;
                symbolMap.put(atom.getSymbol(), val);
            } else {
                symbolMap.put(atom.getSymbol(), 1);
            }
        }
        for (int i = 0; i < target.getAtomCount(); i++) {
            atom = target.getAtom(i);
            atom.setID(Integer.toString(i));
            if (symbolMap.containsKey(atom.getSymbol())) {
                int val = symbolMap.get(atom.getSymbol()) - 1;
                if (val > 0) {
                    symbolMap.put(atom.getSymbol(), val);
                } else {
                    symbolMap.remove(atom.getSymbol());
                }
            }
        }
//        System.out.println("Map " + map);
        return symbolMap.isEmpty();
    }

    private boolean checkAtomMatrix(IAtom queryAtom, IAtom targetAtom) {
        return atomAdjacencyMatrix[Integer.parseInt(queryAtom.getID())][Integer.parseInt(targetAtom.getID())];
    }

    private boolean checkBondMatrix(IBond queryBond, IBond targetBond) {
        return bondAdjacencyMatrix[Integer.parseInt(queryBond.getID())][Integer.parseInt(targetBond.getID())];
    }
}
