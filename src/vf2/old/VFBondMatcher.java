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
package vf2.old;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;

/**
 * Checks if a bond is matching between query and target molecules.
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class VFBondMatcher implements BondMatcher {

    static final long serialVersionUID = -7861469841127328812L;
    private IBond queryBond = null;
    private boolean shouldMatchBonds;

    /**
     * Bond type flag
     */
    /**
     * Constructor
     */
    public VFBondMatcher() {
        this.queryBond = null;
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
        setBondMatchFlag(shouldMatchBonds);
    }

    /** {@inheritDoc}
     * @param targetBond target bond
     * @return true if bonds match
     */
    @Override
    public boolean matches(IBond targetBond) {
        if (!isBondMatchFlag() || (isBondMatchFlag() && isBondTypeMatch(targetBond))) {
            return true;
        }
        return false;
    }

    /**
     * Return true if a bond is matched between query and target
     * @param targetBond
     * @return
     */
    private boolean isBondTypeMatch(IBond targetBond) {
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
