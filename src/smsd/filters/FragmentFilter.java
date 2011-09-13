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
package smsd.filters;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import smsd.AtomAtomMapping;

/**
 * Filter the results based on fragment size.
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 * @cdk.module smsd
 */
public final class FragmentFilter extends BaseFilter implements IChemicalFilter<Integer> {

    private final List<Integer> fragmentSize;

    public FragmentFilter(IAtomContainer rMol, IAtomContainer pMol) {
        super(rMol, pMol);
        fragmentSize = new ArrayList<Integer>();
    }

    @Override
    public synchronized Integer sortResults(
            Map<Integer, AtomAtomMapping> allFragmentAtomMCS,
            Map<Integer, Integer> fragmentScoreMap) throws CDKException {

        int _minFragmentScore = 9999;
        for (Integer Key : allFragmentAtomMCS.keySet()) {
            AtomAtomMapping mcsAtom = allFragmentAtomMCS.get(Key);
            int FragmentCount = getMappedMoleculeFragmentSize(mcsAtom);
            fragmentScoreMap.put(Key, FragmentCount);
            if (_minFragmentScore > FragmentCount) {
                _minFragmentScore = FragmentCount;
            }
        }

        return _minFragmentScore;
    }

    @Override
    public synchronized List<Integer> getScores() {
        return Collections.unmodifiableList(fragmentSize);
    }

    @Override
    public synchronized void clearScores() {
        fragmentSize.clear();
    }

    @Override
    public synchronized void addScore(int counter, Integer value) {
        fragmentSize.add(counter, value);
    }

    @Override
    public synchronized void fillMap(Map<Integer, Integer> fragmentScoreMap) {
        int Index = 0;
        for (Integer score : fragmentSize) {
            fragmentScoreMap.put(Index, score);
            Index++;
        }
    }

    private synchronized int getMappedMoleculeFragmentSize(AtomAtomMapping mcsAtomSolution) {

        IAtomContainer Educt = DefaultChemObjectBuilder.getInstance().newInstance(IMolecule.class, rMol);
        IAtomContainer product = DefaultChemObjectBuilder.getInstance().newInstance(IMolecule.class, pMol);


        if (mcsAtomSolution != null) {
            for (Map.Entry<IAtom, IAtom> map : mcsAtomSolution.getMappings().entrySet()) {
                IAtom atomE = map.getKey();
                IAtom atomP = map.getValue();
                Educt.removeAtomAndConnectedElectronContainers(atomE);
                product.removeAtomAndConnectedElectronContainers(atomP);
            }
        }
        return getFragmentCount(Educt) + getFragmentCount(product);
    }

    private synchronized int getFragmentCount(IAtomContainer molecule) {
        boolean fragmentFlag = true;
        IAtomContainerSet fragmentMolSet = DefaultChemObjectBuilder.getInstance().newInstance(IMoleculeSet.class);
        int countFrag = 0;
        if (molecule.getAtomCount()
                > 0) {
            fragmentFlag = ConnectivityChecker.isConnected(molecule);
            if (!fragmentFlag) {
                fragmentMolSet.add(ConnectivityChecker.partitionIntoMolecules(molecule));
            } else {
                fragmentMolSet.addAtomContainer(molecule);
            }
            countFrag = fragmentMolSet.getAtomContainerCount();
        }
        return countFrag;
    }
}
