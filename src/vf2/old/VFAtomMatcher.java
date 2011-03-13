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

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;

/**
 * Checks if atom is matching between query and target molecules.
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class VFAtomMatcher implements AtomMatcher {

    static final long serialVersionUID = -7861469841127327812L;
    private String symbol = null;
    private IAtom qAtom = null;
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
    public VFAtomMatcher() {
        this.qAtom = null;
        symbol = null;
    }

    /**
     * Constructor
     * @param atom query atom
     * @param shouldMatchBonds bond matching flag
     */
    public VFAtomMatcher(IAtom atom, boolean shouldMatchBonds) {
        this();
        this.qAtom = atom;
        this.symbol = atom.getSymbol();
        setBondMatchFlag(shouldMatchBonds);
    }

    /** {@inheritDoc}
     * @param symbol
     */
    public void setSymbol(String symbol) {
        this.symbol = symbol;
    }

    /** {@inheritDoc}
     *
     * @param targetAtom
     * @return
     */
    @Override
    public boolean matches(IAtom targetAtom) {
        return matchSymbol(targetAtom);
    }

    private boolean matchSymbol(IAtom targetAtom) {
        if (qAtom instanceof IQueryAtom) {
            return ((IQueryAtom) qAtom).matches(targetAtom) ? true : false;
        } else {
            return qAtom.getSymbol().equals(targetAtom.getSymbol()) ? true : false;
        }
    }
}
