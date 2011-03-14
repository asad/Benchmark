/*
 * MX Cheminformatics Tools for Java
 *
 * Copyright (c) 2007-2009 Metamolecular, LLC
 *
 * http://metamolecular.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
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
 */
package isomorphism.vf2.bond;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * This class finds MCS between query and target molecules
 * using VF2 algorithm.
 *
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class VFBondMapper implements IBondMapper {

//    private void setIDs(IAtomContainer atomContainer) {
//        for (int i = 0; i < atomContainer.getAtomCount(); i++) {
//            atomContainer.getAtom(i).setID(String.valueOf(i));
//        }
//    }
    private IAtomContainer query;
    private List<Map<IBond, IBond>> maps;

    /**
     *
     * @param query
     */
    public VFBondMapper(IAtomContainer query) {
//        setIDs(query);
        this.query = query;
        this.maps = new ArrayList<Map<IBond, IBond>>();
    }

    /**
     *
     * @param queryMolecule
     * @param bondMatcher
     */
    public VFBondMapper(IAtomContainer queryMolecule, boolean bondMatcher) {
//        setIDs(queryMolecule);
        this.query = queryMolecule;
        this.maps = new ArrayList<Map<IBond, IBond>>();
    }

    /** {@inheritDoc}
     * @param targetMolecule targetMolecule graph
     */
    @Override
    public boolean hasMap(IAtomContainer targetMolecule) {
//        setIDs(targetMolecule);
        IBondState state = new VFBondState(query, targetMolecule);
        maps.clear();
        boolean flag = mapFirst(state);
        return flag;
    }

    /** {@inheritDoc}
     */
    @Override
    public List<Map<IBond, IBond>> getMaps(IAtomContainer target) {
//        setIDs(target);
        IBondState state = new VFBondState(query, target);
        maps.clear();
        mapAll(state);
        return new ArrayList<Map<IBond, IBond>>(maps);
    }

    /** {@inheritDoc}
     *
     * @param target
     *
     */
    @Override
    public Map<IBond, IBond> getFirstMap(IAtomContainer target) {
//        setIDs(target);
        IBondState state = new VFBondState(query, target);
        mapFirst(state);
        return maps.isEmpty() ? new HashMap<IBond, IBond>() : maps.get(0);
    }

    /** {@inheritDoc}
     */
    @Override
    public int countMaps(IAtomContainer target) {
//        setIDs(target);
        IBondState state = new VFBondState(query, target);
        maps.clear();
        mapAll(state);
        return maps.size();
    }

    private void mapAll(IBondState state) {
        if (state.isDead()) {
            return;
        }

        if (hasMap(state.getMap())) {
            state.backTrack();
        }

        if (state.isGoal()) {
            Map<IBond, IBond> map = state.getMap();
            if (!hasMap(map)) {
                maps.add(state.getMap());
            } else {
                state.backTrack();
            }
        }

        while (state.hasNextCandidate()) {
            VFBondMatcher candidate = state.nextCandidate();
            if (state.isMatchFeasible(candidate)) {
                IBondState nextState = state.nextState(candidate);
                mapAll(nextState);
                nextState.backTrack();
            }
        }
    }

    private boolean mapFirst(IBondState state) {
        if (state.isDead()) {
            return false;
        }

        if (state.isGoal()) {
            maps.add(state.getMap());
            return true;
        }

        boolean found = false;
        while (!found) {
            if (!state.hasNextCandidate()) {
                return false;
            }
            VFBondMatcher candidate = state.nextCandidate();
            if (state.isMatchFeasible(candidate)) {
                IBondState nextState = state.nextState(candidate);
                found = mapFirst(nextState);
                if (found) {
                    return true;
                }
                nextState.backTrack();
            }
        }
        return found;
    }

    private boolean hasMap(Map<IBond, IBond> map) {
        for (Map<IBond, IBond> storedMap : maps) {
            if (storedMap.equals(map)) {
                return true;
            }
        }
        return false;
    }
}