package mcsbenchmark;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import mcsbenchmark.VF2.AtomMapping;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;

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
        List<IMolecule> targets = new ArrayList<IMolecule>();
        IIteratingChemObjectReader tFileReader = read(qFile);
        while (tFileReader.hasNext()) {
            targets.add((IMolecule) tFileReader.next());
        }

        int counter = 0;
        while (qFileReader.hasNext()) {
            IMolecule query = (IMolecule) qFileReader.next();
            int smsdSolutionCount = 0;
            int uitSolutionCount = 0;
            long t0 = System.currentTimeMillis();
            for (IMolecule target : targets) {
                smsdSolutionCount += getSMSDSolutionCount(query, target);
            }
            long timeNow = System.currentTimeMillis();
            long smsdTime = (timeNow - t0);

//            tFileReader = read(tFile);
            long tUIT0 = System.currentTimeMillis();
            for (IMolecule target : targets) {
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

    private static int getSMSDSolutionCount(IMolecule queryMol, IMolecule target) throws CDKException {

//        Substructure substructure = new Substructure();
//        substructure.set(query, target);
//
//        if (substructure.isSubgraph(true)) {
//            return 1;
//        } else {
//            return 0;
//        }

        if (queryMol.getAtomCount() <= target.getAtomCount()) {
            VF2 matcher = new VF2();
            AtomMapping mapping = matcher.isomorphism(queryMol, target);
//            System.out.println("mapping " + mapping);
            if (mapping != null) {
                return 1;
            } else {
                return 0;
            }
        }

        return 0;
    }

    private static int getUITSolutionCount(IMolecule query, IMolecule target) throws CDKException {
//       List bondMapping = UniversalIsomorphismTester.getSubgraphMaps(target, query);
//       List<List<RMap>> sol = UniversalIsomorphismTester.makeAtomsMapsOfBondsMaps(bondMapping, target, query);
        if (UniversalIsomorphismTester.isSubgraph(target, query)) {
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
}
