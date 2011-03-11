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
package vf2.matcher;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;

/**
 * Checks if a bond is matching between query and target molecules.
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.algorithm.vflib.VFLibTest")
public class VFBondMatcher implements BondMatcher {

    static final long serialVersionUID = -7861469841127328812L;
    private IBond queryBond = null;
    private int unsaturation = 0;
    private boolean shouldMatchBonds;
    private IQueryBond smartQueryBond = null;

    /**
     * Bond type flag
     */
    /**
     * Constructor
     */
    public VFBondMatcher() {
        this.queryBond = null;
        this.unsaturation = -1;
        shouldMatchBonds = false;
    }

    /**
     * Constructor
     * @param queryMol query Molecule
     * @param queryBond query Molecule
     * @param shouldMatchBonds bond match flag
     */
    public VFBondMatcher(IAtomContainer queryMol, IBond queryBond, boolean shouldMatchBonds) {
        super();
        this.queryBond = queryBond;
//        this.unsaturation = getUnsaturation(queryMol, this.queryBond);
        setBondMatchFlag(shouldMatchBonds);
    }

    /**
     * Constructor
     * @param queryBond query Molecule
     */
    public VFBondMatcher(IQueryBond queryBond) {
        super();
        this.smartQueryBond = queryBond;
    }

    /** {@inheritDoc}
     *
     * @param targetConatiner target container
     * @param targetBond target bond
     * @return true if bonds match
     */
    @Override
    public boolean matches(IAtomContainer targetConatiner, IBond targetBond) {
        if (this.smartQueryBond != null) {
            return smartQueryBond.matches(targetBond);
        } else if (!isBondMatchFlag() || (isBondMatchFlag() && isBondTypeMatch(targetBond))) {
            return true;
        }
//        if (this.unsaturation != -1 && this.unsaturation == getUnsaturation(targetConatiner, targetBond)) {
//            return true;
//        }
        return false;
    }

    /**
     * Return true if a bond is matched between query and target
     * @param targetBond
     * @return
     */
    private boolean isBondTypeMatch(IBond targetBond) {
        if ((queryBond.getFlag(CDKConstants.ISAROMATIC) == targetBond.getFlag(CDKConstants.ISAROMATIC))
                && (queryBond.getOrder() == targetBond.getOrder())) {
            return true;
        } else if (queryBond.getFlag(CDKConstants.ISAROMATIC) && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
            return true;
        }
        return false;
    }

    private int getUnsaturation(IAtomContainer container, IBond bond) {
        return getUnsaturation(container, bond.getAtom(0)) + getUnsaturation(container, bond.getAtom(1));
    }

    private int getUnsaturation(IAtomContainer container, IAtom atom) {
        return getValency(atom) - container.getConnectedAtomsCount(atom);
    }

    private int getValency(IAtom atom) {
        return (atom.getValency() == null) ? 0 : atom.getValency().intValue();
    }

//    private int getUnsaturation(IAtomContainer container, IAtom atom) {
//        return getValency(atom) - (countNeighbors(container, atom) + countImplicitHydrogens(atom));
//    }
//    private int countNeighbors(IAtomContainer container, IAtom atom) {
//        return container.getConnectedAtomsCount(atom);
//    }
//
//    private int countImplicitHydrogens(IAtom atom) {
//        return (atom.getImplicitHydrogenCount() == null)
//                ? 0 : atom.getImplicitHydrogenCount();
//    }
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
}
