package test;

import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;

import vf2.AtomMapping;
import vf2.VF2;

public class PermutationTest {
    
    private void addAtoms(IMolecule mol, int count, String symbol) {
        for (int i = 0; i < count; i++) {
            mol.addAtom(new Atom(symbol));
        }
    }
    
    @Test
    public void cyclopentane() throws CDKException {
        IMolecule cycA = new Molecule();
        addAtoms(cycA, 5, "C");
        cycA.addBond(0, 1, IBond.Order.SINGLE);
        cycA.addBond(0, 2, IBond.Order.SINGLE);
        cycA.addBond(1, 4, IBond.Order.SINGLE);
        cycA.addBond(2, 3, IBond.Order.SINGLE);
        cycA.addBond(3, 4, IBond.Order.SINGLE);
        
        IMolecule cycB = new Molecule();
        addAtoms(cycB, 5, "C");
        cycB.addBond(0, 3, IBond.Order.SINGLE);
        cycB.addBond(0, 4, IBond.Order.SINGLE);
        cycB.addBond(1, 2, IBond.Order.SINGLE);
        cycB.addBond(1, 3, IBond.Order.SINGLE);
        cycB.addBond(2, 4, IBond.Order.SINGLE);
        
        VF2 vf2 = new VF2();
        AtomMapping mapping = vf2.isomorphism(cycA, cycB);
        System.out.println(mapping);
        
        boolean uitThinks = UniversalIsomorphismTester.isIsomorph(cycA, cycB);
        System.out.println("UIT thinks " + uitThinks);
    }

}
