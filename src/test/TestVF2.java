package test;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;

public class TestVF2 {

    public static void main(String[] args) throws InvalidSmilesException, CDKException {
        TestVF2Coverage vfTest = new TestVF2Coverage();
        vfTest.TestCovergare();
        PermutationTest pt = new PermutationTest();
        pt.cyclopentane();
    }
}
