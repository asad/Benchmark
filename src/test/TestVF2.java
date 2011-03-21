/*
 * TestIsomorphism.java
 *
 * Created on January 28, 2007, 2:06 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package test;

import chemkit.AtomMapping;
import chemkit.VF2;
import helper.ImageGenerator;
import helper.GraphMolecule;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.IChemObjectReader;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class TestVF2 {

    /**
     * 
     * @param molFile
     * @return
     */
    public static IMolecule readMol(String molFile) {
        IMolecule mol = null;
        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        try {
            MDLReader reader = new MDLReader(new FileReader(molFile), IChemObjectReader.Mode.RELAXED);
            IMolecule molObject = (IMolecule) reader.read(new GraphMolecule());
            reader.close();
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molObject);
            CDKHueckelAromaticityDetector.detectAromaticity(molObject);
            sdg.setMolecule(new GraphMolecule(molObject));
            sdg.generateCoordinates();
            mol = sdg.getMolecule();
            setID(mol);
        } catch (Exception ex) {
            Logger.getLogger(TestVF2.class.getName()).log(Level.SEVERE, null, ex);
        }
        return mol;
    }

    private static void setID(IMolecule mol) {
        int index = 1;
        for (IAtom atom : mol.atoms()) {
            atom.setID(String.valueOf(index++));
        }
    }
    static private ImageGenerator imageGenerator = null;

    /** Creates a new instance of TestIsomorphism */
    public TestVF2() {
    }

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws Exception 
     */
    public static void main(String[] args) throws IOException, CDKException, Exception {
//        TestVF2Coverage vfTest = new TestVF2Coverage();
//        vfTest.TestCovergare();
//        vfTest.TestCovergareOldVF2();
//        PermutationTest pt = new PermutationTest();
//        pt.cyclopentane();
//        pt.cyclobutane();
//        pt.permuteSubgraph();
//        pt.permuteSubgraph2();
//        pt.permuteSubgraph3();
//        pt.permuteSubgraph4();

        testIsomorphism();
    }

    public static void testIsomorphism() throws Exception {

        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());


//        String target = "CC(C)(COP([O-])([O-])(O*)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C1N=CN=C2N)C(O)C(=O)NCCC(=O)NCCS";
//        String query = "CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C1N=CN=C2N)C(O)C(=O)NCCC(=O)NCCS";
////////
//        String query = "Nc1ccccc1";
////        String target = "C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C";
////        //Phosphate matches
//        String query = "[H]OC1C(COP([O-])(=O)OP([O-])(*)=O)OC(C1O)n1cnc2c(N)ncnc12";
//        String target = "[H]C(O)OC1C(COP([O-])(=O)OP([O-])(*)=O)OC(C1O)n1cnc2c(N)ncnc12";
//
////        String query= "[O-][CH]=[O]";//"CC(=O)C=O";//"O=C(CN1C=NC=CC1=O)NCC2=CC=CC=C2";
////        String target= "[O-][CH]=[O]";//"C[C@@H](O)C-O";//"CC(C)(C)N1C2=C(C=N1)C(=O)N(C=N2)CC(=O)NCC3=CC=CC=C3Cl";
//
//        //String query = "[H][O][C@]([H])([C]([H])=[O])[C]([H])([H])[H]";
//        String query = "[H][O][C]([H])([H])[C](=[O])[C]([H])([H])[O][P]([O-])([O-])=[O]";
//        String target = "[H][O][C@]1([H])[C@]([H])([O][P]([O-])([O-])=[O])[O][C@@]([H])([C]([H])([H])[H])[C@@]([H])([O][H])[C@@]1([H])[O][H]";
//
//        String query="[H]C([H])(C([O-])=O)C([H])([H])C(=O)SCCN";
//        String target ="[H]C([H])([H])[C@@]([H])(C([O-])=O)C(=O)SCC";
//
//        String target = "[H]C([H])([H])[C@@]([H])(C([O-])=O)C(=O)S[*]";
//        String query = "[H]C([H])(C([O-])=O)C([H])([H])C(=O)S[*]";
////
////        String query = //"Nc1ncnc2n(cnc12)[C@H]1O[C@@H](COP([O-])(=O)OP([O-])=O)[C@H](OP([O-])([O-])=O)[C@@H]1O";
////                "[H]OP(O)([O-])([O-])OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12";
//////                "Nc1ncnc2n(cnc12)[C@H]1O[C@@H](COP([O-])(=O)OP(O)([O-])=O)[C@H](OP([O-])([O-])=O)[C@@H]1O";
////        String target = "[H]OP(O)([O-])([O-])OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12";
//
//
//        String query = "[H][O][C+]([O][P]([O-])(=[O])[O][C]([H])([H])[C]1([H])[O][C]([H])([n]2[c]([H])[n][c]3[c]([n][c]"
//                + "([H])[n][c]23)[N]([H])[H])[C]([H])([O][H])[C]1([H])[O][H])[c]1[c]([H])[c]([H])[c]([H])[n+]([c]1[H])"
//                + "[C]1([H])[O][C]([H])([C]([H])([H])[O][P]([O-])(=[O])[O][P]([O-])(=[O])[O][C]([H])([H])[C]2([H])[O][C]"
//                + "([H])([n]3[c]([H])[n][c]4[c]([n][c]([H])[n][c]34)[N]([H])[H])[C@]([H])([O][H])[C@]2([H])[O][H])[C@@]"
//                + "([H])([O][H])[C@@]1([H])[O][H]";
//        String target = "[H][O][C@@]1([H])[C]([H])([O][C]([H])([C]([H])([H])[O][P]([O-])(=[O])[O][P]([O-])(=[O])[O][C]"
//                + "([H])([H])[C]2([H])[O][C]([H])([n+]3[c]([H])[c]([H])[c]([H])[c]([c]3[H])[C](=[O])[N]([H])[H])[C@]"
//                + "([H])([O][H])[C@]2([H])[O][H])[C@@]1([H])[O][H])[n]1[c]([H])[n][c]2[c]([n][c]([H])[n][c]12)[N]([H])[H]";
//////
//
////        //This is from MACiE reaction 192.02
//        String query = "[H][N](*)[C](=[O])[C]([H])([N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]1([H])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](*)=[O])[C]([H])([H])[S][S]*)[C]([H])([H])[S][S][C]1([H])[H])[C]([H])([H])[S-]";
//        String target = "[H][N](*)[C](=[O])[C]1([H])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](*)=[O])[C]([H])([H])[S][S]*)[C]([H])([H])[S-])[C]([H])([H])[S][S][C]1([H])[H]";

//        /*Ring matches*/
//        String target = "[H][O][C@@]([H])([C]([H])([H])[H])[C@]([H])([O][H])[C@]1([H])[N]([H])[C]2([O][H])[C]([O-])=[N][C](=[N][C]2=[N][C]1([H])[H])[N]([H])[H]";
//        String query = "[H][O][C@@]([H])([C]([H])([H])[H])[C@]([H])([O][H])[C@]1([H])[N+]([H])=[c]2[c]([O-])[n][c]([n][c]2=[N][C]1([H])[H])[N]([H])[H]";
////

//        String target = "C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C";
//        String query = "Nc1ccccc1";

//        //Subgraph match
//
        String query = "CC";
        String target = "C1CCC12CCCC2";
//        //Single MApping test
////        String target = "O";
////        String query = "Nc1ncnc2n(cnc12)C1OC(COP([O-])=O)C(O)C1O.Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])(=O)OCC2OC([C@H](O)[C@@H]2O)[n+]2cccc([CH+]O)c2)[C@@H](O)[C@H]1O";
//
//
//        String query = "C1CC1";
//        String target = "CC(C)C";
//
//        String query = "CC";
////        String query="C1CC1";//"C1CC1(CC1CC1)";//Gillian
//        String query = "C1CCCCC1";//Hexane
//         String query = "C1=CC=CC=C1";//BENZENE
////        String target = "C1CCC12CCCC2";
//        String target = "C1CCCCC1";//Hexane//
//        String target = "C1=CC=CC=C1";//BENZENE
////        String target = "C1CC1(CC1CC1)";//"C1CC1";//gillian

//        String query = "CCCc1nn(C)c2c1nc(nc2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1";//cdk bug expected 9 got 8
//        String target = "NC1=NC2=C(N=CN2[C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)C(=O)N1";//cdk bug expected 9 got 8
////
//        String query = "O=C(CN1C=NC=CC1=O)NCC2=CC=CC=C2";//rguha query
//        String target = "CC(C)(C)N1C2=C(C=N1)C(=O)N(C=N2)CC(=O)NCC3=CC=CC=C3Cl";

//        //Sameer hassan
//        String query = "(c1cc(c[n+](c1)[C@H]2[C@@H]([C@@H]([C@H](O2)CO[P@@](=O)([O-])O[P@@](=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)n4cnc5c4ncnc5N)O)O)O)O)C(=O)N)";//sameer query
//        String target = "(C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C4N=CN=C5N)O)O)O)O)";



        IAtomContainer mol1 = new GraphMolecule(sp.parseSmiles(query));
//        int counter = 0;
//        System.out.println("Atoms: ");
//        for (IAtom atom : mol1.atoms()) {
//            atom.setID(String.valueOf(counter++));
//            System.out.println(atom.getSymbol() + " " + atom.getID());
//        }
//        counter = 0;
//        for (IBond bond : mol1.bonds()) {
//            bond.setID(String.valueOf(counter++));
//        }
//        System.out.println("Cloned Atoms: ");
//        IAtomContainer molClone = (AtomContainer) mol1.clone();
//        for (IAtom atom : molClone.atoms()) {
//            System.out.println(atom.getSymbol() + " " + atom.getID());
//        }
        IAtomContainer mol2 = new GraphMolecule(sp.parseSmiles(target));
        printMolecules(mol1);
        printMolecules(mol2);

        VF2 comparison = new VF2();
        List<AtomMapping> maps = comparison.isomorphisms(mol1, mol2, true);

        System.out.println("Mol1 Size. " + mol1.getAtomCount());
        System.out.println("Mol2 Size. " + mol2.getAtomCount());
        System.out.println("Map Size. " + maps.size());

        printMapping(maps, mol1, mol2);
        generateImage("Subgraph", mol1, mol2, maps);
    }

    private static void printMapping(List<AtomMapping> comparison, IAtomContainer query, IAtomContainer target) {

        int count_final_sol = 0;
        System.out.println("Output of the final Mappings: ");
        try {
            if (!comparison.isEmpty()) {

                for (AtomMapping final_solution : comparison) {
                    int final_solution_size = final_solution.getSize();
                    System.out.println("Final mapping Nr. " + ++count_final_sol + " Size:" + final_solution_size);

                    for (Map.Entry<IAtom, IAtom> mapping : final_solution.getAtomMapping().entrySet()) {
                        System.out.println(query.getAtomNumber(mapping.getKey()) + " " + target.getAtomNumber(mapping.getValue()));

                        IAtom eAtom = mapping.getKey();
                        IAtom pAtom = mapping.getValue();

                        System.out.println(eAtom.getSymbol() + "[" + eAtom.getFormalCharge() + "] "
                                + pAtom.getSymbol() + "[" + pAtom.getFormalCharge() + "] ");
                    }
                    System.out.println("");
                }

                System.out.println("");
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private static void generateImage(String outPutFileName, IAtomContainer query, IAtomContainer target, List<AtomMapping> smsd) throws Exception {

        imageGenerator = new ImageGenerator();

        ////set the format right for the Tanimoto score (only two digits printed)
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(2);
        nf.setMinimumFractionDigits(2);
        System.out.println("Output of the final Mappings: ");
        int counter = 1;
        for (AtomMapping mapping : smsd) {
            String tanimoto = "0";
            String stereo = "NA";
            String label = "Scores [" + "Tanimoto: " + tanimoto + ", Stereo: " + stereo + "]";
            Map<Integer, Integer> map = new HashMap<Integer, Integer>(mapping.getSize());
            for (Map.Entry<IAtom, IAtom> mapAtom : mapping.getAtomMapping().entrySet()) {
//                System.out.println(query.getAtomNumber(mapAtom.getKey()) + " " + target.getAtomNumber(mapAtom.getValue()));
                map.put(query.getAtomNumber(mapAtom.getKey()), target.getAtomNumber(mapAtom.getValue()));
            }

            imageGenerator.addImages(query, target, label, map);
            counter++;
        }
        String filePNG = System.getProperty("user.dir") + File.separator + outPutFileName;
        imageGenerator.createImage(filePNG, "Query", "Target");
    }

    /**
     * 
     * @param Molecule1
     * @param Molecule2
     */
    private static void printMolecules(IAtomContainer Molecule1, IAtomContainer Molecule2) {

        System.out.println("Molecule 1");

        for (int i = 0; i < Molecule1.getAtomCount(); i++) {

            System.out.print(Molecule1.getAtom(i).getSymbol() + " ");
        }

        System.out.println();
        System.out.println("Molecule 2");
        for (int i = 0; i < Molecule2.getAtomCount(); i++) {

            System.out.print(Molecule2.getAtom(i).getSymbol() + " ");
        }
        System.out.println();

    }

    private static void printMolecules(IAtomContainer Molecule1) {

        System.out.println("Molecule 1: " + Molecule1.getAtomCount());

        for (int i = 0; i < Molecule1.getAtomCount(); i++) {

            System.out.print(Molecule1.getAtom(i).getSymbol());
//            System.out.print(Molecule1.getAtom(i).getSymbol() + "[" + Molecule1.getAtom(i).getProperty(CANONICAL_LABEL) + "]");
        }

        System.out.println();
    }
}
