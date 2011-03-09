package mcsbenchmark;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import vf2.AtomMapping;
import vf2.VF2;

public class TestVF2 {

    public static void main(String[] args) throws InvalidSmilesException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CCOC(=O)c1[nH]c2ccc(C)cc2(c1(N))");
        IAtomContainer target = sp.parseSmiles("CCOC(=O)c2[nH]c1ccc(C)cc1c2(N=CN(CC)CC)");
////        IAtomContainer target = sp.parseSmiles("CCOC(=O)c1[nH]c2ccc(C)cc2(c1(N))");
//        IAtomContainer query = sp.parseSmiles("CCO");
//        IAtomContainer target = sp.parseSmiles("CCOCN");
        if (query.getAtomCount() <= target.getAtomCount()) {
            VF2 matcher = new VF2();
            AtomMapping mapping = matcher.isomorphism(query, target);
            System.out.println("mapping " + mapping);
        }
    }
}
