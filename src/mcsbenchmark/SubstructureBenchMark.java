package mcsbenchmark;

import helper.CDKSMILES;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import mcsbenchmark.VF2.AtomMapping;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Asad
 */
public class SubstructureBenchMark {

    /**
     * @param args the command line arguments
     * @throws FileNotFoundException
     * @throws IOException
     * @throws CDKException  
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, CDKException {
        String queryFilePath = (args.length > 0) ? args[0] : "data/actives.sdf";
        String targetFilePath = (args.length > 1) ? args[1] : "data/all.sdf";
        File qFile = new File(queryFilePath);
        File tFile = new File(targetFilePath);

        IIteratingChemObjectReader qFileReader = read(qFile);

        if (qFileReader == null) {
            throw new IOException("Unknown input type ");
        }
        List<IAtomContainer> targets = new ArrayList<IAtomContainer>();
        IIteratingChemObjectReader tFileReader = read(tFile);
        while (tFileReader.hasNext()) {
            IAtomContainer ac = configure((IMolecule) tFileReader.next());
            targets.add(ac);
        }

        int counter = 0;
        while (qFileReader.hasNext()) {
            IAtomContainer query = (IMolecule) qFileReader.next();
            query = configure(query);
            int smsdSolutionCount = 0;
            int uitSolutionCount = 0;
            long t0 = System.currentTimeMillis();
            for (IAtomContainer target : targets) {
                smsdSolutionCount += getSMSDSolutionCount(query, target);
            }
            long timeNow = System.currentTimeMillis();
            long smsdTime = (timeNow - t0);

//            tFileReader = read(tFile);
            long tUIT0 = System.currentTimeMillis();
            for (IAtomContainer target : targets) {
                uitSolutionCount += getUITSolutionCount(query, target);
            }
            timeNow = System.currentTimeMillis();
            long uitTime = timeNow - tUIT0;

            String out = String.format("%d SMSDt %d SMSDs %d UITt %d UITs %d ",
                    counter, smsdTime, smsdSolutionCount, uitTime, uitSolutionCount);
            System.out.println(out);
            counter++;
        }

    }

    private static int getSMSDSolutionCount(IAtomContainer query, IAtomContainer target) throws CDKException {

//        Substructure substructure = new Substructure();
//        substructure.set(query, target);
//
//        if (substructure.isSubgraph(true)) {
//            return 1;
//        } else {
//            return 0;
//        }

        if (query.getAtomCount() <= target.getAtomCount()) {
            VF2 matcher = new VF2();
//            CDKSMILES cdkSmiles = new CDKSMILES();
//            String a = cdkSmiles.getSMILES(query);
//            String b = cdkSmiles.getSMILES(target);
//            AtomMapping mapping = matcher.isomorphism(a, b);
//            System.out.println("mapping " + mapping);
            AtomMapping mapping = matcher.isomorphism(query, target);
            if (!mapping.isEmpty()) {
//                
//                System.out.println("\nVF");
//                System.out.println(cdkSmiles.getSMILES(query));
//                System.out.println(cdkSmiles.getSMILES(target));
//                System.out.println();
                return 1;
            } else {
                return 0;
            }
        }

        return 0;
    }

    private static int getUITSolutionCount(IAtomContainer query, IAtomContainer target) throws CDKException {
//       List bondMapping = UniversalIsomorphismTester.getSubgraphMaps(target, query);
//       List<List<RMap>> sol = UniversalIsomorphismTester.makeAtomsMapsOfBondsMaps(bondMapping, target, query);
        if (UniversalIsomorphismTester.isSubgraph(target, query)) {
//            CDKSMILES cdkSmiles = new CDKSMILES();
//            System.out.println("\nUIT");
//            System.out.println(cdkSmiles.getSMILES(query));
//            System.out.println(cdkSmiles.getSMILES(target));
//            System.out.println();
            return 1;
        } else {
            return 0;
        }
//       return sol.size();
    }

    /**
     *
     * @param file
     * @return
     * @throws FileNotFoundException
     */
    public static IIteratingChemObjectReader read(File file) throws FileNotFoundException {

        FileReader in = new FileReader(file);
        return new IteratingMDLReader(
                in, NoNotificationChemObjectBuilder.getInstance());

    }

    private static IAtomContainer configure(IAtomContainer atomContainer) throws CDKException {
//        /**
//         * Prepare the target molecule for analysis. <p/> We perform ring perception and aromaticity detection and set up
//         * the appropriate properties. Right now, this function is called each time we need to do a query and this is
//         * inefficient.
//         *
//         * @throws CDKException if there is a problem in ring perception or aromaticity detection, which is usually related
//         *                      to a timeout in the ring finding code.
//         */
//        // Code copied from
//        // org.openscience.cdk.qsar.descriptors.atomic.AtomValenceDescriptor;
//        Map<String, Integer> valencesTable = new HashMap<String, Integer>();
//        valencesTable.put("H", 1);
//        valencesTable.put("Li", 1);
//        valencesTable.put("Be", 2);
//        valencesTable.put("B", 3);
//        valencesTable.put("C", 4);
//        valencesTable.put("N", 5);
//        valencesTable.put("O", 6);
//        valencesTable.put("F", 7);
//        valencesTable.put("Na", 1);
//        valencesTable.put("Mg", 2);
//        valencesTable.put("Al", 3);
//        valencesTable.put("Si", 4);
//        valencesTable.put("P", 5);
//        valencesTable.put("S", 6);
//        valencesTable.put("Cl", 7);
//        valencesTable.put("K", 1);
//        valencesTable.put("Ca", 2);
//        valencesTable.put("Ga", 3);
//        valencesTable.put("Ge", 4);
//        valencesTable.put("As", 5);
//        valencesTable.put("Se", 6);
//        valencesTable.put("Br", 7);
//        valencesTable.put("Rb", 1);
//        valencesTable.put("Sr", 2);
//        valencesTable.put("In", 3);
//        valencesTable.put("Sn", 4);
//        valencesTable.put("Sb", 5);
//        valencesTable.put("Te", 6);
//        valencesTable.put("I", 7);
//        valencesTable.put("Cs", 1);
//        valencesTable.put("Ba", 2);
//        valencesTable.put("Tl", 3);
//        valencesTable.put("Pb", 4);
//        valencesTable.put("Bi", 5);
//        valencesTable.put("Po", 6);
//        valencesTable.put("At", 7);
//        valencesTable.put("Fr", 1);
//        valencesTable.put("Ra", 2);
//        valencesTable.put("Cu", 2);
//        valencesTable.put("Mn", 2);
//        valencesTable.put("Co", 2);
//
//        // do all ring perception
//        AllRingsFinder arf = new AllRingsFinder();
//        IRingSet allRings;
//        try {
//            allRings = arf.findAllRings(atomContainer);
//        } catch (CDKException e) {
////            logger.debug(e.toString());
//            throw new CDKException(e.toString(), e);
//        }
//
//        // sets SSSR information
//        SSSRFinder finder = new SSSRFinder(atomContainer);
//        IRingSet sssr = finder.findEssentialRings();
//
//        for (IAtom atom : atomContainer.atoms()) {
//
//            // add a property to each ring atom that will be an array of
//            // Integers, indicating what size ring the given atom belongs to
//            // Add SSSR ring counts
//            if (allRings.contains(atom)) { // it's in a ring
//                atom.setFlag(CDKConstants.ISINRING, true);
//                // lets find which ring sets it is a part of
//                List<Integer> ringsizes = new ArrayList<Integer>();
//                IRingSet currentRings = allRings.getRings(atom);
//                int min = 0;
//                for (int i = 0; i < currentRings.getAtomContainerCount(); i++) {
//                    int size = currentRings.getAtomContainer(i).getAtomCount();
//                    if (min > size) {
//                        min = size;
//                    }
//                    ringsizes.add(size);
//                }
//                atom.setProperty(CDKConstants.RING_SIZES, ringsizes);
//                atom.setProperty(CDKConstants.SMALLEST_RINGS, sssr.getRings(atom));
//            } else {
//                atom.setFlag(CDKConstants.ISINRING, false);
//            }
//
//            // determine how many rings bonds each atom is a part of
//            int hCount;
//            if (atom.getImplicitHydrogenCount() == CDKConstants.UNSET) {
//                hCount = 0;
//            } else {
//                hCount = atom.getImplicitHydrogenCount();
//            }
//
//            List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);
//            int total = hCount + connectedAtoms.size();
//            for (IAtom connectedAtom : connectedAtoms) {
//                if (connectedAtom.getSymbol().equals("H")) {
//                    hCount++;
//                }
//            }
//            atom.setProperty(CDKConstants.TOTAL_CONNECTIONS, total);
//            atom.setProperty(CDKConstants.TOTAL_H_COUNT, hCount);
//
//            if (valencesTable.get(atom.getSymbol()) != null) {
//                int formalCharge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();
//                atom.setValency(valencesTable.get(atom.getSymbol()) - formalCharge);
//            }
//        }
//
//        for (IBond bond : atomContainer.bonds()) {
//            if (allRings.getRings(bond).getAtomContainerCount() > 0) {
//                bond.setFlag(CDKConstants.ISINRING, true);
//            }
//        }
//
//        for (IAtom atom : atomContainer.atoms()) {
//            List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);
//
//            int counter = 0;
//            IAtom any;
//            for (IAtom connectedAtom : connectedAtoms) {
//                any = connectedAtom;
//                if (any.getFlag(CDKConstants.ISINRING)) {
//                    counter++;
//                }
//            }
//            atom.setProperty(CDKConstants.RING_CONNECTIONS, counter);
//        }
//
//        // check for atomaticity
//        try {
//            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
//            CDKHueckelAromaticityDetector.detectAromaticity(atomContainer);
//        } catch (CDKException e) {
////            logger.debug(e.toString());
//            throw new CDKException(e.toString(), e);
//        }

        return atomContainer;
    }
}
