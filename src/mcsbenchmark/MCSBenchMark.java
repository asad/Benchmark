package mcsbenchmark;

import com.metamolecular.mx.io.daylight.SMILESReader;
import com.metamolecular.mx.map.DefaultMapper;
import com.metamolecular.mx.map.Mapper;
import com.metamolecular.mx.model.Molecule;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;

/**
 *
 * @author Asad
 */
public class MCSBenchMark {

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
        List<Molecule> targetsMX = new ArrayList<Molecule>();
        IIteratingChemObjectReader tFileReader = read(tFile);
        while (tFileReader.hasNext()) {
            targets.add((IMolecule) tFileReader.next());
        }

        int counter = 0;
        while (qFileReader.hasNext()) {
            IMolecule query = (IMolecule) qFileReader.next();
            int smsdSolutionCount = 0;
            int uitSolutionCount = 0;
            int mxSolutionCount = 0;
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

            Molecule queryMX = getMXMolecule(query);
            for (IMolecule target : targets) {
                Molecule targetMX = getMXMolecule(target);
                targetsMX.add(targetMX);
            }

            long tMX0 = System.currentTimeMillis();
            for (Molecule target : targetsMX) {
                mxSolutionCount += getMXSolutionCount(queryMX, target);
            }
            timeNow = System.currentTimeMillis();
            long mxTime = timeNow - tMX0;



            String out = String.format("%d SMSDt %d SMSDs %d UITt %d UITs %d MXt %d MXs %d",
                    counter, smsdTime, smsdSolutionCount, uitTime, uitSolutionCount,
                    mxTime, mxSolutionCount);
            System.out.println(out);
            counter++;
        }

    }

    private static int getSMSDSolutionCount(IMolecule query, IMolecule target) throws CDKException {

        SubStructure substructure = new SubStructure();
        substructure.set(query, target);

        if (substructure.isSubgraph(true)) {
            return 1;
        } else {
            return 0;
        }
    }

    private static int getMXSolutionCount(Molecule query, Molecule target) throws CDKException {

        Mapper mapper = new DefaultMapper(query);
        if (mapper.hasMap(target)) {
            return 1;
        } else {
            return 0;
        }
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

    public static Molecule getMXMolecule(IAtomContainer cdkMolecule) throws CDKException, IOException {
        SmilesGenerator sg = new SmilesGenerator(true);
        String smiles = "";
        smiles = sg.createSMILESWithoutCheckForMultipleMolecules(cdkMolecule, false, new boolean[cdkMolecule.getBondCount()]);
        System.out.println("SMILE " + smiles);
        Molecule mxMol = SMILESReader.read(smiles);
        return mxMol;
    }
}
