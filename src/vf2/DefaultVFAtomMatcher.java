/*
 * Copyright (C) 2006-2010  Syed Asad Rahman <asad@ebi.ac.uk>
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
package vf2;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;

/**
 * Checks if atom is matching between query and target molecules.
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.algorithm.vflib.VFLibTest")
public class DefaultVFAtomMatcher implements VFAtomMatcher {

    static final long serialVersionUID = -7861469841127327812L;
    private int maximumNeighbors;
    private String symbol = null;
    private IAtom qAtom = null;
    private IQueryAtom smartQueryAtom = null;
    private boolean shouldMatchBonds = false;

    /**
     * @return the shouldMatchBonds
     */
    public boolean isBondMatchFlag() {
        return shouldMatchBonds;
    }

    /**
     * @param shouldMatchBonds the shouldMatchBonds to set
     */
    public final void setBondMatchFlag(boolean shouldMatchBonds) {
        this.shouldMatchBonds = shouldMatchBonds;
    }

    /**
     * Constructor
     */
    public DefaultVFAtomMatcher() {
        this.qAtom = null;
        symbol = null;
        maximumNeighbors = -1;
    }

    /**
     * Constructor
     * @param queryContainer query atom container
     * @param atom query atom
     * @param shouldMatchBonds bond matching flag
     */
    public DefaultVFAtomMatcher(IAtomContainer queryContainer, IAtom atom, boolean shouldMatchBonds) {
        this();
        this.qAtom = atom;
        this.symbol = atom.getSymbol();
        setBondMatchFlag(shouldMatchBonds);
//        this.maximumNeighbors = countImplicitHydrogens(atom)
//                + queryContainer.getConnectedAtomsCount(atom);

//        System.out.println("Atom " + atom.getSymbol());
//        System.out.println("MAX allowed " + maximumNeighbors);
    }

    /**
     * Constructor
     * @param smartQueryAtom query atom
     * @param container 
     */
    public DefaultVFAtomMatcher(IQueryAtom smartQueryAtom, IQueryAtomContainer container) {
        this();
        this.smartQueryAtom = smartQueryAtom;
        this.symbol = smartQueryAtom.getSymbol();
    }

    /**
     * Constructor
     * @param queryContainer query atom container
     * @param template query atom
     * @param blockedPositions
     * @param shouldMatchBonds bond matching flag
     */
    public DefaultVFAtomMatcher(IAtomContainer queryContainer, IAtom template, int blockedPositions, boolean shouldMatchBonds) {
        this(queryContainer, template, shouldMatchBonds);
        this.maximumNeighbors = countImplicitHydrogens(template)
                + queryContainer.getConnectedAtomsCount(template)
                - blockedPositions;
    }

    /**
     *
     * @param maximum numbers of connected atoms allowed
     */
    public void setMaximumNeighbors(int maximum) {
        this.maximumNeighbors = maximum;
    }

    /** {@inheritDoc}
     * @param symbol
     */
    public void setSymbol(String symbol) {
        this.symbol = symbol;
    }

    /** {@inheritDoc}
     *
     * @param targetContainer
     * @param targetAtom
     * @return
     */
    @Override
    public boolean matches(IAtomContainer targetContainer, IAtom targetAtom) {
        if (smartQueryAtom != null && qAtom == null) {
            if (!smartQueryAtom.matches(targetAtom)) {
                return false;
            }
        } else if (!matchSymbol(targetAtom)) {
            return false;
        }

//        if (!matchMaximumNeighbors(targetContainer, targetAtom)) {
//            return false;
//        }
        return true;
    }

    private boolean matchSymbol(IAtom atom) {
        if (symbol == null) {
            return false;
        }
        return symbol.equals(atom.getSymbol());
    }

    private boolean matchMaximumNeighbors(IAtomContainer targetContainer, IAtom targetAtom) {
        if (maximumNeighbors == -1 || !isBondMatchFlag()) {
            return true;
        }

        int maximumTargetNeighbors = targetContainer.getConnectedAtomsCount(targetAtom)
                + countImplicitHydrogens(targetAtom);
        return maximumTargetNeighbors <= maximumNeighbors;
    }

    private int countImplicitHydrogens(IAtom atom) {
        return (atom.getImplicitHydrogenCount() == null)
                ? 0 : atom.getImplicitHydrogenCount();
    }
}
