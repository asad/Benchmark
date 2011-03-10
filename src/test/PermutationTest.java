package test;

import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.smiles.SmilesParser;
import tools.labelling.AtomContainerAtomPermutor;
import tools.labelling.AtomContainerPrinter;

import vf2.AtomMapping;
import vf2.VF2;

public class PermutationTest {

    private void addAtoms(IMolecule mol, int count, String symbol) {
        for (int i = 0; i < count; i++) {
            mol.addAtom(new Atom(symbol));
        }
    }

    private void testWithUIT(IAtomContainer molA, IAtomContainer molB) throws CDKException {
        boolean uitThinks = UniversalIsomorphismTester.isIsomorph(molA, molB);
        System.out.println("UIT thinks " + uitThinks);
    }

    private void testWithVF2(IAtomContainer molA, IAtomContainer molB) {
        VF2 vf2 = new VF2();
        AtomMapping mapping = vf2.isomorphism(molA, molB);
        System.out.println(mapping);
    }

    @Test
    public void cyclobutane() throws CDKException {
        IMolecule cycA = new Molecule();
        addAtoms(cycA, 4, "C");
        cycA.addBond(0, 1, IBond.Order.SINGLE);
        cycA.addBond(0, 2, IBond.Order.SINGLE);
        cycA.addBond(1, 3, IBond.Order.SINGLE);
        cycA.addBond(2, 3, IBond.Order.SINGLE);

        IMolecule cycB = new Molecule();
        addAtoms(cycB, 4, "C");
        cycB.addBond(0, 2, IBond.Order.SINGLE);
        cycB.addBond(0, 3, IBond.Order.SINGLE);
        cycB.addBond(1, 2, IBond.Order.SINGLE);
        cycB.addBond(1, 3, IBond.Order.SINGLE);

        testWithVF2(cycA, cycB);
        testWithUIT(cycA, cycB);
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

        testWithVF2(cycA, cycB);
        testWithUIT(cycA, cycB);
    }

    @Test
    public void permuteSubgraph() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer q = sp.parseSmiles("CCOC(=O)c1[nH]c2ccc(C)cc2(c1(N))");
        IAtomContainer t = sp.parseSmiles("CCOC(=O)c2[nH]c1ccc(C)cc1c2(N=CN(CC)CC)");

        AtomContainerPrinter printer = new AtomContainerPrinter();


        System.out.println("\n Step 1");

        System.out.println("Input Q: " + printer.toString(q));
        System.out.println("Input T: " + printer.toString(t));

        testWithVF2(q, t);
        testWithUIT(t, q);

        AtomContainerAtomPermutor acpQuery = new AtomContainerAtomPermutor(q);
        AtomContainerAtomPermutor acpTarget = new AtomContainerAtomPermutor(t);

        System.out.println("\n Step 2");


        IAtomContainer query = acpQuery.next();
        IAtomContainer target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 3");


        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 4");

        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 5");


        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 6");

        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 7");

        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

    }

    @Test
    public void permuteSubgraph2() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer q = sp.parseSmiles("NC1=CC=C(O)[N+]2=C1C=CC=C2");
        IAtomContainer t = sp.parseSmiles("NC1=CC(P(O)C2=CNC=C2)=C(O)[N+]2=C1C=CC=C2");

        AtomContainerPrinter printer = new AtomContainerPrinter();


        System.out.println("\n Step 1");

        System.out.println("Input Q: " + printer.toString(q));
        System.out.println("Input T: " + printer.toString(t));

        testWithVF2(q, t);
        testWithUIT(t, q);

        AtomContainerAtomPermutor acpQuery = new AtomContainerAtomPermutor(q);
        AtomContainerAtomPermutor acpTarget = new AtomContainerAtomPermutor(t);

        System.out.println("\n Step 2");


        IAtomContainer query = acpQuery.next();
        IAtomContainer target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 3");


        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 4");

        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 5");


        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 6");

        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

        System.out.println("\n Step 7");

        query = acpQuery.next();
        target = acpTarget.next();

        System.out.println("Input Q: " + printer.toString(query));
        System.out.println("Input T: " + printer.toString(target));

        testWithVF2(query, target);
        testWithUIT(target, query);

    }
}
