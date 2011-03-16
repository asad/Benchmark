package mcsbenchmark;

import helper.Molecule;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import loop.SMSDVF2;
import loop.UITLoop;
import loop.ChemkitVF2;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;

/**
 *
 * @author Asad
 */
public class SubstructureBenchMark {

    /**
     * @param args the command line arguments
     * @throws FileNotFoundException
     * @throws Exception  
     */
    public static void main(String[] args) throws FileNotFoundException, Exception {
        String queryFilePath = (args.length > 0) ? args[0] : "data/actives.sdf";
        String targetFilePath = (args.length > 1) ? args[1] : "data/all.sdf";
//        
//        String queryFilePath = (args.length > 0) ? args[0] : "data/q.sdf";
//        String targetFilePath = (args.length > 1) ? args[1] : "data/t.sdf";

//        String queryFilePath = (args.length > 0) ? args[0] : "data/1Query.sdf";
//        String targetFilePath = (args.length > 1) ? args[1] : "data/4Targets.sdf";

//        String queryFilePath = (args.length > 0) ? args[0] : "data/some.sdf";
//        String targetFilePath = (args.length > 1) ? args[1] : "data/some.sdf";



        File qFile = new File(queryFilePath);
        File tFile = new File(targetFilePath);

        IIteratingChemObjectReader qFileReader = read(qFile);

        if (qFileReader == null) {
            throw new IOException("Unknown input type ");
        }
        List<IMolecule> queries = new ArrayList<IMolecule>();
        while (qFileReader.hasNext()) {
            IMolecule mol = new Molecule((IMolecule) qFileReader.next());
            queries.add(mol);
        }
//        Collections.shuffle(queries);

        List<IMolecule> targets = new ArrayList<IMolecule>();
        IIteratingChemObjectReader tFileReader = read(tFile);

        int counter = 1;
        while (tFileReader.hasNext()) {
            IMolecule mol = new Molecule((IMolecule) tFileReader.next());
            mol.setID(String.valueOf(counter++));
            targets.add(mol);
        }
        counter = 1;
        System.out.println(
                "number of queries " + queries.size()
                + " number of targets " + targets.size());

        for (IMolecule query : queries) {
            String out = String.format("%d ", counter);
            out += String.format("\t%d ", query.getAtomCount());

            query.setID(String.valueOf(counter));

            UITLoop uitLoop = new UITLoop();
            uitLoop.run(query, targets);
            out += uitLoop;

            SMSDVF2 smsdLoop = new SMSDVF2();
            smsdLoop.run(query, targets);
            out += smsdLoop;

            ChemkitVF2 vfLibLoop = new ChemkitVF2();
            vfLibLoop.run(query, targets);
            out += vfLibLoop;

//            MCSPlusLoop mcsPlusLoop = new MCSPlusLoop();
//            mcsPlusLoop.run(query, targets);
//            out += mcsPlusLoop;

            System.out.println(out);
            counter++;
        }

    }

    /**
     *
     * @param file
     * @return
     * @throws FileNotFoundException
     */
    public static IIteratingChemObjectReader read(File file) throws Exception {
        FileReader in = new FileReader(file);
        String path = file.getAbsolutePath();
        if (path.endsWith(".sdf")) {
            return new IteratingMDLReader(
                    in, NoNotificationChemObjectBuilder.getInstance());
        } else if (path.endsWith(".smi")) {
            return new IteratingSMILESReader(
                    in, NoNotificationChemObjectBuilder.getInstance());
        } else {
            throw new Exception("Unrecognised filetype " + path);
        }

    }

    private static IAtomContainer configure(IAtomContainer atomContainer) throws CDKException {
        return atomContainer;
    }
}
