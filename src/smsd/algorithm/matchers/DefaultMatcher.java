/* Copyright (C) 2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
package smsd.algorithm.matchers;

import org.openscience.cdk.interfaces.IBond;

/**
 * Checks if atom is matching between query and target molecules.
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class DefaultMatcher {

    /**
     * 
     * @param bondMatcher
     * @param bondA2
     * @param shouldMatchBonds
     * @return
     */
    private static boolean isBondMatch(BondMatcher queryBondMatcher, IBond targetBond) {
        return queryBondMatcher.matches(targetBond);
    }

    private static boolean isAtomMatch(IBond bondA1, IBond bondA2, boolean shouldMatchRings) {

        AtomMatcher atomMatcher1 = new DefaultAtomMatcher(bondA1.getAtom(0), shouldMatchRings);
        AtomMatcher atomMatcher2 = new DefaultAtomMatcher(bondA1.getAtom(1), shouldMatchRings);

        // ok, atoms match
        if (atomMatcher1.matches(bondA2.getAtom(0)) && atomMatcher2.matches(bondA2.getAtom(1))) {
//            System.out.println("Atom Matched");
            return true;
        }
        // ok, atoms match
        if (atomMatcher1.matches(bondA2.getAtom(1)) && atomMatcher2.matches(bondA2.getAtom(0))) {
//            System.out.println("Atom Matched");
            return true;
        }
        return false;
    }

    /**
     * 
     * @param bondA1
     * @param bondA2
     * @param matchBond
     * @param shouldMatchRings
     * @return
     */
    public static boolean matches(IBond bondA1, IBond bondA2, boolean matchBond, boolean shouldMatchRings) {

        if (!isAtomMatch(bondA1, bondA2, shouldMatchRings)) {
            return false;
        }
        if (matchBond && !isBondMatch(new DefaultBondMatcher(bondA1, matchBond), bondA2)) {
            return false;
        }
        return true;
    }
}
