package mcsbenchmark;

import gui.ImageGenerator;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.smsd.algorithm.mcsplus.MCSPlusHandler;
import org.openscience.smsd.algorithm.rgraph.CDKMCSHandler;
import org.openscience.smsd.algorithm.vflib.VFlibMCSHandler;
import org.openscience.smsd.tools.MolHandler;

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
    public static void main(String[] args) throws FileNotFoundException, IOException, CDKException, Exception {
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
        int counter = 1;

        while (qFileReader.hasNext()) {
            Map<Integer, List<Integer>> mapSol = new HashMap<Integer, List<Integer>>();
            IMolecule query = (IMolecule) qFileReader.next();

            int index = 1;
            for (IMolecule target : targets) {
                List<Integer> l = new ArrayList<Integer>(6);
                l.add(0, 0);
                l.add(1, 0);
                l.add(2, 0);
                l.add(3, 0);
                l.add(4, 0);
                l.add(5, 0);
                mapSol.put(index, l);
                MDLV2000Writer writer = new MDLV2000Writer(
                        new FileWriter(new File(Integer.toString(index++) + ".mol")));
                writer.write(target);
                writer.close();
            }


            int smsdSolutionCount = 0;
            int uitSolutionCount = 0;
            long t0 = System.currentTimeMillis();
            int tCounter = 1;
            for (IMolecule target : targets) {
//                System.out.println("ID " + target.getID());
                smsdSolutionCount += getSMSDSolutionCount(tCounter, query, target, mapSol);
                tCounter++;
            }
//            
            long timeNow = System.currentTimeMillis();
            long smsdTime = (timeNow - t0);

//            tFileReader = read(tFile);
            long tUIT0 = System.currentTimeMillis();
            tCounter = 1;
            for (IMolecule target : targets) {
                uitSolutionCount += getUITSolutionCount(tCounter, query, target, mapSol);
                tCounter++;
            }
            timeNow = System.currentTimeMillis();
            long uitTime = timeNow - tUIT0;
//
//            long tMCSPlus0 = System.currentTimeMillis();
//            tCounter = 1;
//            for (IMolecule target : targets) {
//                uitSolutionCount += getMCSPlusSolutionCount(tCounter, query, target, mapSol);
//                tCounter++;
//            }
//            timeNow = System.currentTimeMillis();
//            long MCSPlusTime = timeNow - tMCSPlus0;

            String out = String.format("%d SMSDt %d SMSDs %d UITt %d UITs %d ",
                    counter, smsdTime, smsdSolutionCount, uitTime, uitSolutionCount);
            System.out.println(out);
            System.out.println();
            int counterT = 1;
            for (List<Integer> l : mapSol.values()) {
                String outT = String.format("%d SMSD (Sol Count) %d SMSD (Sol Size) %d UIT (Sol Count) %d UIT (Sol Size) %d ",
                        counterT++, l.get(0), l.get(1), l.get(2), l.get(3));
                System.out.println(outT);
            }

            counter++;
            break;
        }

    }

    private static int getSMSDSolutionCount(int counter, IMolecule query, IMolecule target, Map<Integer, List<Integer>> mapSol) throws CDKException, Exception {

        VFlibMCSHandler mcs = null;
        mcs = new VFlibMCSHandler();
        mcs.set(new MolHandler(query, false, false), new MolHandler(target, false, false));
        mcs.searchMCS(true);
        if (mcs.getFirstMapping() != null && !mcs.getFirstMapping().isEmpty()) {
            List<Map<Integer, Integer>> map = mcs.getAllMapping();
            if (mapSol.containsKey(counter)) {
                int solcount = map.size();
                int solsize = mcs.getFirstAtomMapping().size();
                List<Integer> l = mapSol.get(counter);
                l.set(0, solcount);
                l.set(1, solsize);
                mapSol.put(counter, l);
            }
            generateImage(counter + "_SMSD", query, target, map);
            return 1;
        } else {
            return 0;
        }
    }

    private static int getUITSolutionCount(int counter, IMolecule query, IMolecule target, Map<Integer, List<Integer>> mapSol) throws CDKException, Exception {
//       List bondMapping = UniversalIsomorphismTester.getSubgraphMaps(target, query);
//       List<List<RMap>> sol = UniversalIsomorphismTester.makeAtomsMapsOfBondsMaps(bondMapping, target, query);
//        List<IAtomContainer> list = UniversalIsomorphismTester.getOverlaps(query, target);

        CDKMCSHandler cdkmcs = new CDKMCSHandler();
        cdkmcs.set(new MolHandler(query, false, false), new MolHandler(target, false, false));
        cdkmcs.searchMCS(true);

        if (cdkmcs.getFirstMapping() != null && !cdkmcs.getFirstMapping().isEmpty()) {
            List<Map<Integer, Integer>> map = cdkmcs.getAllMapping();
            if (mapSol.containsKey(counter)) {
                int solcount = map.size();
                int solsize = cdkmcs.getFirstAtomMapping().size();
                List<Integer> l = mapSol.get(counter);
                l.set(2, solcount);
                l.set(3, solsize);
                mapSol.put(counter, l);
            }
            generateImage(counter + "_UIT", query, target, map);
            return 1;
        } else {
            return 0;
        }
//       return sol.size();
    }

    private static int getMCSPlusSolutionCount(int counter, IMolecule query, IMolecule target, Map<Integer, List<Integer>> mapSol) throws CDKException, Exception {
//       List bondMapping = UniversalIsomorphismTester.getSubgraphMaps(target, query);
//       List<List<RMap>> sol = UniversalIsomorphismTester.makeAtomsMapsOfBondsMaps(bondMapping, target, query);
//        List<IAtomContainer> list = UniversalIsomorphismTester.getOverlaps(query, target);

        MCSPlusHandler mcsplus = new MCSPlusHandler();
        mcsplus.set(new MolHandler(query, false, false), new MolHandler(target, false, false));
        mcsplus.searchMCS(true);

        if (mcsplus.getFirstMapping() != null && !mcsplus.getFirstMapping().isEmpty()) {
            List<Map<Integer, Integer>> map = mcsplus.getAllMapping();
            if (mapSol.containsKey(counter)) {
                int solcount = map.size();
                int solsize = mcsplus.getFirstAtomMapping().size();
                List<Integer> l = mapSol.get(counter);
                l.set(4, solcount);
                l.set(5, solsize);
                mapSol.put(counter, l);
            }
            generateImage(counter + "_MCSPlus", query, target, map);
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

    private static void generateImage(String outPutFileName, IAtomContainer query, IAtomContainer target, List<Map<Integer, Integer>> smsd) throws Exception {

        ImageGenerator imageGenerator = new ImageGenerator();

        ////set the format right for the Tanimoto score (only two digits printed)
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(2);
        nf.setMinimumFractionDigits(2);
        System.out.println("Output of the final Mappings: ");
        int counter = 1;
        for (Map<Integer, Integer> mapping : smsd) {
            String tanimoto = mapping.size() + "";
            String stereo = "NA";
            String label = "Scores [" + "Size: " + tanimoto + ", Stereo: " + stereo + "]";
            imageGenerator.addImages(query, target, label, mapping);
            counter++;
        }
        String filePNG = System.getProperty("user.dir") + File.separator + outPutFileName;
        imageGenerator.createImage(filePNG, "Query", "Target");

    }
}
