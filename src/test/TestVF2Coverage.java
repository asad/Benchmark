/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import vf2.AtomMapping;
import vf2.VF2;
import vf2.old.VFMapper;

/**
 *
 * @author Asad
 */
public class TestVF2Coverage {

    String[] queries = {"c1ccccc1C(N)C(O)=O",
        "CC(C)OCC(C)=C",
        "Nc1ccc(CC)cc1",
        "O=c1[nH]c(=O)c2ccccc12",
        "O=C4C=C3CCC2C1CCCC1CCC2C3CC4",
        "CC2(C)CN1C(=O)C(N)C1S2",
        "OC1CCCCOCC1"};

    private List<IAtomContainer> readSMILES() {
        List<IAtomContainer> list = new ArrayList<IAtomContainer>();

        return list;
    }

    @Test
    public void TestCovergare() throws InvalidSmilesException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CCOC(=O)c1[nH]c2ccc(C)cc2(c1(N))");
        IAtomContainer target = sp.parseSmiles("CCOC(=O)c2[nH]c1ccc(C)cc1c2(N=CN(CC)CC)");
        if (query.getAtomCount() <= target.getAtomCount()) {
            VF2 matcher = new VF2();
            AtomMapping mapping = matcher.isomorphism(query, target);
            System.out.println("mapping " + mapping);

//            List<AtomMapping> mapping = matcher.isomorphisms(query, target);
//            for (AtomMapping map : mapping) {
//                System.out.println("mapping " + map);
//            }
        }
    }

    @Test
    public void TestCovergareOldVF2() throws InvalidSmilesException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CCOC(=O)c1[nH]c2ccc(C)cc2(c1(N))");
        IAtomContainer target = sp.parseSmiles("CCOC(=O)c2[nH]c1ccc(C)cc1c2(N=CN(CC)CC)");
        if (query.getAtomCount() <= target.getAtomCount()) {
            VFMapper matcher = new VFMapper(query);
            if (matcher.hasMap(target)) {
                System.out.println("mapping " + matcher.getFirstMap(target).size());
            }
        }
    }
}
