package test;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;

public class TestVF2 {

    /**
     * 
     * @param args
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    public static void main(String[] args) throws InvalidSmilesException, CDKException {
        TestVF2Coverage vfTest = new TestVF2Coverage();
        vfTest.TestCovergare();
        vfTest.TestCovergareOldVF2();
        PermutationTest pt = new PermutationTest();
        pt.cyclopentane();
        pt.cyclobutane();
        pt.permuteSubgraph();
    }
}
