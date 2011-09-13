/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import chemkit.substructure.VF2;
import java.util.ArrayList;
import java.util.List;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import smsd.algorithm.vflib.VF2Sub;


/**
 *
 * @author Asad
 */
public class VF2CoverageTest {

    @Test
    public void TestCovergare() throws InvalidSmilesException, CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CCOC(=O)c1[nH]c2ccc(C)cc2(c1(N))");
        IAtomContainer target = sp.parseSmiles("CCOC(=O)c2[nH]c1ccc(C)cc1c2(N=CN(CC)CC)");
        if (query.getAtomCount() <= target.getAtomCount()) {
            VF2 matcher = new VF2(true, true);
            matcher.set(query, target);
            if (matcher.isSubgraph()) {
                System.out.println("mapping " + matcher);
            }
        }
    }

    @Test
    public void TestCovergareOldVF2() throws InvalidSmilesException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CCOC(=O)c1[nH]c2ccc(C)cc2(c1(N))");
        IAtomContainer target = sp.parseSmiles("CCOC(=O)c2[nH]c1ccc(C)cc1c2(N=CN(CC)CC)");
        if (query.getAtomCount() <= target.getAtomCount()) {
            VF2Sub matcher = new VF2Sub(true, true);
            matcher.set(query, target);
            if (matcher.isSubgraph()) {
                System.out.println("mapping " + matcher);
            }
        }
    }
    
    
//    @Before
//    String[] queries = {"c1ccccc1C(N)C(O)=O",
//        "CC(C)OCC(C)=C",
//        "Nc1ccc(CC)cc1",
//        "O=c1[nH]c(=O)c2ccccc12",
//        "O=C4C=C3CCC2C1CCCC1CCC2C3CC4",
//        "CC2(C)CN1C(=O)C(N)C1S2",
//        "OC1CCCCOCC1"};
//
//    private List<IAtomContainer> readSMILES() {
//        List<IAtomContainer> list = new ArrayList<IAtomContainer>();
//        return list;
//    }
}
