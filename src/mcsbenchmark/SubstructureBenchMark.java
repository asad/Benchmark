package mcsbenchmark;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import vf2.AtomMapping;
import vf2.VF2;

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
        String queryFilePath = (args.length > 0) ? args[0] : "data/some.sdf";//"data/actives.sdf";
        String targetFilePath = (args.length > 1) ? args[1] : "data/some.sdf";//"data/t.sdf";

//        String queryFilePath = (args.length > 0) ? args[0] : "data/q.sdf";
//        String targetFilePath = (args.length > 1) ? args[1] : "data/t.sdf";

        File qFile = new File(queryFilePath);
        File tFile = new File(targetFilePath);

        IIteratingChemObjectReader qFileReader = read(qFile);

        if (qFileReader == null) {
            throw new IOException("Unknown input type ");
        }
        List<IAtomContainer> targets = new ArrayList<IAtomContainer>();
        IIteratingChemObjectReader tFileReader = read(tFile);
        while (tFileReader.hasNext()) {
            IAtomContainer ac = configure((IAtomContainer) tFileReader.next());
            targets.add(ac);
        }

        int counter = 0;
        while (qFileReader.hasNext()) {
            IAtomContainer query = (IAtomContainer) qFileReader.next();
            query = configure(query);

//            AtomContainerPrinter printer = new AtomContainerPrinter();
//            System.out.println(printer.toString(query));


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
        return atomContainer;
    }
}
