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
package vf2.old;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

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

    /**
     * initialize the VFState with query and target
     * @param query
     * @param target
     */
    public VFState(IAtomContainer query, IAtomContainer target) {
        this.map = new AtomMapping(target, query);
        this.queryPath = new ArrayList<IAtom>();
        this.targetPath = new ArrayList<IAtom>();

        this.query = query;
        this.target = target;
        this.candidates = new ArrayList<Match>();
        loadRootCandidates();
    }

    private VFState(VFState state, Match match) {
        this.candidates = new ArrayList<Match>();
        this.queryPath = new ArrayList<IAtom>(state.queryPath);
        this.targetPath = new ArrayList<IAtom>(state.targetPath);

        this.map = state.map;
        this.query = state.query;
        this.target = state.target;

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
        if (!matchAtoms(match)) {
            return false;
        }
        if (!matchBonds(match)) {
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

    private void loadRootCandidates() {
        for (int i = 0; i < query.getAtomCount(); i++) {
            IAtom qAtom = query.getAtom(i);
            for (int j = 0; j < target.getAtomCount(); j++) {
                IAtom tAtom = target.getAtom(j);
                Match match = new Match(qAtom, tAtom);
                if (matchAtoms(match)) {
                    candidates.add(match);
                }
            }
        }
        System.out.println("Compatibility graph " + candidates.size());
    }

//@TODO Asad Check the Neighbour count
    private void loadCandidates(Match lastMatch) {
        IAtom atom = lastMatch.getTargetAtom();
        List<IAtom> targetNeighbors = target.getConnectedAtomsList(atom);
        List<IAtom> queryNeighbors = query.getConnectedAtomsList(lastMatch.getQueryAtom());
        for (IAtom queryAtom : queryNeighbors) {
            for (IAtom targetAtom : targetNeighbors) {
                Match match = new Match(queryAtom, targetAtom);
                if (matchAtoms(match)) {
                    if (candidateFeasible(match)) {
                        candidates.add(match);
                    }
                }
            }
        }
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

    private boolean matchAtoms(Match match) {
        IAtom targetAtom = match.getTargetAtom();
        if (query.getConnectedAtomsCount(match.getQueryAtom()) > target.getConnectedAtomsCount(targetAtom)) {
            return false;
        }
        return matchAtoms(match.getQueryAtom(), targetAtom);
    }

    private boolean matchBonds(Match match) {
        if (queryPath.isEmpty()) {
            return true;
        }

        if (!matchBondsToHead(match)) {
            return false;
        }

        for (int i = 0; i < queryPath.size() - 1; i++) {
            IBond queryBond = query.getBond(queryPath.get(i), match.getQueryAtom());
            IBond targetBond = target.getBond(targetPath.get(i), match.getTargetAtom());
            if (queryBond == null) {
                continue;
            }

            if (targetBond == null) {
                return false;
            }
            if (!matchBond(queryBond, targetBond)) {
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
        return matchBond(queryBond, targetBond);
    }

    private IAtom getQueryPathHead() {
        return queryPath.get(queryPath.size() - 1);
    }

    private IAtom getTargetPathHead() {
        return targetPath.get(targetPath.size() - 1);
    }

    boolean matchBond(IBond sourceBond, IBond targetBond) {
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

    private boolean matchAtomTypes1(IBond qbond, IBond tbond) {
        if (qbond.getAtom(0).getSymbol().equals(tbond.getAtom(0).getSymbol())
                && qbond.getAtom(1).getSymbol().equals(tbond.getAtom(1).getSymbol())) {
            return true;
        }
        return false;
    }

    private boolean matchAtomTypes2(IBond qbond, IBond tbond) {
        if (qbond.getAtom(0).getSymbol().equals(tbond.getAtom(1).getSymbol())
                && qbond.getAtom(1).getSymbol().equals(tbond.getAtom(0).getSymbol())) {
            return true;
        }
        return false;
    }
}
