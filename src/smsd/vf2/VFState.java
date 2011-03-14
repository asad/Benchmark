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
package smsd.vf2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;

/**
 * This class finds mapping states between query and target
 * molecules.
 * 
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class VFState implements IState {

    private List<Match> candidates;
    private IAtomContainer query;
    private IAtomContainer target;
    private List<IAtom> queryPath;
    private List<IAtom> targetPath;
    private AtomMapping map;
    private Map<IAtom, List<IAtom>> neighbourQueryMap;
    private Map<IAtom, List<IAtom>> neighbourTargetMap;
    private boolean[][] atomAdjacencyMatrix;
    private boolean[][] bondAdjacencyMatrix;

    /**
     * initialize the VFState with query and target
     * @param query
     * @param target
     */
    public VFState(IAtomContainer query, IAtomContainer target) {
        this.query = query;
        this.target = target;
        this.map = new AtomMapping(target, query);
        this.candidates = new ArrayList<Match>();


        this.neighbourQueryMap = new HashMap<IAtom, List<IAtom>>();
        this.neighbourTargetMap = new HashMap<IAtom, List<IAtom>>();
//        initialize(query, target);

        if (testIsSubgraphHeuristics(query, target)) {
            this.queryPath = new ArrayList<IAtom>(query.getAtomCount());
            this.targetPath = new ArrayList<IAtom>(target.getAtomCount());
            this.atomAdjacencyMatrix = new boolean[query.getAtomCount()][target.getAtomCount()];
            this.bondAdjacencyMatrix = new boolean[query.getBondCount()][target.getBondCount()];
            if (loadRootCandidates()) {
                for (IBond bondQ : query.bonds()) {
                    for (IBond bondT : target.bonds()) {
                        if (matchBond(bondQ, bondT)) {
                            bondAdjacencyMatrix[Integer.parseInt(bondQ.getID())][Integer.parseInt(bondT.getID())] = true;
                        } else {
                            bondAdjacencyMatrix[Integer.parseInt(bondQ.getID())][Integer.parseInt(bondT.getID())] = false;
                        }
                    }
                }
            }
        }
    }

    private VFState(VFState state, Match match) {
        this.candidates = new ArrayList<Match>();
        this.queryPath = new ArrayList<IAtom>(state.queryPath);
        this.targetPath = new ArrayList<IAtom>(state.targetPath);

        this.map = state.map;
        this.query = state.query;
        this.target = state.target;

        this.neighbourQueryMap = state.neighbourQueryMap;
        this.neighbourTargetMap = state.neighbourTargetMap;
        this.atomAdjacencyMatrix = state.atomAdjacencyMatrix;
        this.bondAdjacencyMatrix = state.bondAdjacencyMatrix;

        map.add(match.getQueryAtom(), match.getTargetAtom());
        queryPath.add(match.getQueryAtom());
        targetPath.add(match.getTargetAtom());
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
    public Map<IAtom, IAtom> getMap() {
        return new HashMap<IAtom, IAtom>(map.getMapping());
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
        return query.getAtomCount() > target.getAtomCount();
    }

    /** {@inheritDoc}
     */
    @Override
    public boolean isGoal() {
        return map.size() == query.getAtomCount();
    }

    /** {@inheritDoc}
     */
    @Override
    public boolean isMatchFeasible(Match match) {
        if (map.containsQueryAtom(match.getQueryAtom())
                || map.containsTargetAtom(match.getTargetAtom())) {
            return false;
        }
//        if (!atomMatcher(match)) {
//            return false;
//        }
        if (!checkAtomMatrix(match.getQueryAtom(), match.getTargetAtom())) {
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
    public Match nextCandidate() {
        return candidates.remove(candidates.size() - 1);
    }

    /** {@inheritDoc}
     */
    @Override
    public IState nextState(Match match) {
        return new VFState(this, match);
    }

    private boolean loadRootCandidates() {
        for (IAtom qAtom : neighbourQueryMap.keySet()) {
            boolean flag = false;
            for (IAtom tAtom : neighbourTargetMap.keySet()) {
                Match match = new Match(qAtom, tAtom);
                if (atomMatcher(match)) {
                    candidates.add(match);
                    atomAdjacencyMatrix[Integer.parseInt(qAtom.getID())][Integer.parseInt(tAtom.getID())] = true;
                    flag = true;
                } else {
                    atomAdjacencyMatrix[Integer.parseInt(qAtom.getID())][Integer.parseInt(tAtom.getID())] = false;
                }
            }
            if (!flag) {
                candidates.clear();
                return false;
            }
        }
        return true;
//        System.out.println("Compatibility graph " + candidates.size());
    }

//@TODO Asad Check the Neighbour count
    private void loadCandidates(Match lastMatch) {
        List<IAtom> queryNeighbors = neighbourQueryMap.get(lastMatch.getQueryAtom());
        List<IAtom> targetNeighbors = neighbourTargetMap.get(lastMatch.getTargetAtom());

        for (IAtom queryAtom : queryNeighbors) {
            for (IAtom targetAtom : targetNeighbors) {
                if (checkAtomMatrix(queryAtom, targetAtom)) {
                    Match match = new Match(queryAtom, targetAtom);
                    if (candidateFeasible(match)) {
//                    System.out.println("map " + map.size());
                        candidates.add(match);
                    }
                }
            }
        }
//        System.out.println("candidates " + candidates.size());
    }

    private boolean candidateFeasible(Match candidate) {
        for (IAtom queryAtom : map.queryAtoms()) {
            if (queryAtom.equals(candidate.getQueryAtom())
                    || map.getMappedTargetAtom(queryAtom).equals(candidate.getTargetAtom())) {
                return false;
            }
        }
        return true;
    }
    //This function is updated by Asad to include more matches

    private boolean atomMatcher(Match match) {
        if (neighbourQueryMap.get(match.getQueryAtom()).size() > neighbourTargetMap.get(match.getTargetAtom()).size()) {
            return false;
        }
        return matchAtoms(match.getQueryAtom(), match.getTargetAtom());
    }

    private boolean bondMatcher(Match match) {
        if (queryPath.isEmpty()) {
            return true;
        }

        if (!matchBondsToHead(match)) {
            return false;
        }

        for (int i = 0; i < queryPath.size() - 1; i++) {
            IBond queryBond = query.getBond(queryPath.get(i), match.getQueryAtom());
            if (queryBond == null) {
                continue;
            }
            IBond targetBond = target.getBond(targetPath.get(i), match.getTargetAtom());
            if (targetBond == null) {
                return false;
            }
//            if (!matchBond(queryBond, targetBond)) {
//                return false;
//            }
            if (!checkBondMatrix(queryBond, targetBond)) {
                return false;
            }
        }
        return true;
    }

    private boolean isHeadMapped() {
        IAtom head = queryPath.get(queryPath.size() - 1);
        List<IAtom> queryHeadNeighbors = query.getConnectedAtomsList(head);
        for (IAtom neighbor : queryHeadNeighbors) {
            if (!map.containsQueryAtom(neighbor)) {
                return false;
            }
        }
        return true;
    }

    private boolean matchBondsToHead(Match match) {
        IAtom queryHead = getQueryPathHead();
        IAtom targetHead = getTargetPathHead();

        IBond queryBond = query.getBond(queryHead, match.getQueryAtom());
        IBond targetBond = target.getBond(targetHead, match.getTargetAtom());

        if (queryBond == null || targetBond == null) {
            return false;
        }
//        return matchBond(queryBond, targetBond);
        return checkBondMatrix(queryBond, targetBond);
    }

    private IAtom getQueryPathHead() {
        return queryPath.get(queryPath.size() - 1);
    }

    private IAtom getTargetPathHead() {
        return targetPath.get(targetPath.size() - 1);
    }

    private boolean matchBond(IBond queryBond, IBond targetBond) {
        if (queryBond instanceof IQueryBond) {
            return ((IQueryBond) queryBond).matches(targetBond);
        } else if ((queryBond.getFlag(CDKConstants.ISAROMATIC) == targetBond.getFlag(CDKConstants.ISAROMATIC))
                && (queryBond.getOrder() == targetBond.getOrder())) {
            return true;
        } else if (queryBond.getFlag(CDKConstants.ISAROMATIC) && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
            return true;
        }
        return false;
    }

    boolean matchAtoms(IAtom sourceAtom, IAtom targetAtom) {
        if (sourceAtom instanceof IQueryAtom) {
            return ((IQueryAtom) sourceAtom).matches(targetAtom) ? true : false;
        } else {
            return sourceAtom.getSymbol().equals(targetAtom.getSymbol()) ? true : false;
        }
    }

    /**
     *  Checks some simple heuristics for whether the subgraph query can
     *  realistically be atom subgraph of the supergraph. If, for example, the
     *  number of nitrogen atoms in the query is larger than that of the supergraph
     *  it cannot be part of it.
     *
     * @param  ac1  the subgraph to be checked. 
     * @param  ac2  the super-graph to be tested for. Must not be an IQueryAtomContainer.
     * @return    true if the subgraph ac1 has atom chance to be atom subgraph of ac2
     * @throws org.openscience.cdk.exception.CDKException if the first molecule is an instance
     * of IQueryAtomContainer
     */
    private boolean testIsSubgraphHeuristics(IAtomContainer ac1, IAtomContainer ac2) {

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
            bond.setID(Integer.toString(i));
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
            bond.setID(Integer.toString(i));
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
        for (int i = 0; i < ac1.getAtomCount(); i++) {
            atom = ac1.getAtom(i);
            atom.setID(Integer.toString(i));
            neighbourQueryMap.put(atom, query.getConnectedAtomsList(atom));
            if (symbolMap.containsKey(atom.getSymbol())) {
                int val = symbolMap.get(atom.getSymbol()) + 1;
                symbolMap.put(atom.getSymbol(), val);
            } else {
                symbolMap.put(atom.getSymbol(), 1);
            }
        }
        for (int i = 0; i < ac2.getAtomCount(); i++) {
            atom = ac2.getAtom(i);
            atom.setID(Integer.toString(i));
            neighbourTargetMap.put(atom, target.getConnectedAtomsList(atom));
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
//    /**
//     *  Checks some simple heuristics for whether the subgraph query can
//     *  realistically be atom subgraph of the supergraph. If, for example, the
//     *  number of nitrogen atoms in the query is larger than that of the supergraph
//     *  it cannot be part of it.
//     *
//     * @param  ac1  the subgraph to be checked. 
//     * @param  ac2  the super-graph to be tested for. Must not be an IQueryAtomContainer.
//     * @return    true if the subgraph ac1 has atom chance to be atom subgraph of ac2
//     * @throws org.openscience.cdk.exception.CDKException if the first molecule is an instance
//     * of IQueryAtomContainer
//     */
//    private void initialize(IAtomContainer ac1, IAtomContainer ac2) {
//
//        IBond bond = null;
//
//        for (int i = 0; i < ac1.getBondCount(); i++) {
//            bond = ac1.getBond(i);
//            bond.setID(Integer.toString(i));
//        }
//        for (int i = 0; i < ac2.getBondCount(); i++) {
//            bond = ac2.getBond(i);
//            bond.setID(Integer.toString(i));
//        }
//
//        IAtom atom = null;
//        for (int i = 0; i < ac1.getAtomCount(); i++) {
//            atom = ac1.getAtom(i);
//            atom.setID(Integer.toString(i));
//            neighbourQueryMap.put(atom, query.getConnectedAtomsList(atom));
//        }
//        for (int i = 0; i < ac2.getAtomCount(); i++) {
//            atom = ac2.getAtom(i);
//            atom.setID(Integer.toString(i));
//            neighbourTargetMap.put(atom, target.getConnectedAtomsList(atom));
//        }
//    }

    private boolean checkAtomMatrix(IAtom queryAtom, IAtom targetAtom) {
        return atomAdjacencyMatrix[Integer.parseInt(queryAtom.getID())][Integer.parseInt(targetAtom.getID())];
    }

    private boolean checkBondMatrix(IBond queryBond, IBond targetBond) {
        return bondAdjacencyMatrix[Integer.parseInt(queryBond.getID())][Integer.parseInt(targetBond.getID())];
    }
}
