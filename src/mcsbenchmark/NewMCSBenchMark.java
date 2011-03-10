package mcsbenchmark;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import loop.MCSPlusLoop;
import loop.SMSDLoop;
import loop.UITLoop;
import loop.VFLibLoop;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;

/**
 *
 * @author Asad
 */
public class NewMCSBenchMark {

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
        List<IMolecule> queries = new ArrayList<IMolecule>();
        while (qFileReader.hasNext()) {
            queries.add((IMolecule) qFileReader.next());
        }
        Collections.shuffle(queries);
        
        List<IMolecule> targets = new ArrayList<IMolecule>();
        IIteratingChemObjectReader tFileReader = read(tFile);
        while (tFileReader.hasNext()) {
            targets.add((IMolecule) tFileReader.next());
        }
        int counter = 1;
        System.out.println(
                "number of queries " + queries.size() +
                " number of targets " + targets.size());

        for (IMolecule query : queries) {
            String out = String.format("%d ", counter);

            UITLoop uitLoop = new UITLoop();
            uitLoop.run(query, targets);
            out += uitLoop;
            
            SMSDLoop smsdLoop = new SMSDLoop();
            smsdLoop.run(query, targets);
            out += smsdLoop;

            VFLibLoop vfLibLoop = new VFLibLoop();
            vfLibLoop.run(query, targets);
            out += vfLibLoop;

            MCSPlusLoop mcsPlusLoop = new MCSPlusLoop();
            mcsPlusLoop.run(query, targets);
            out += mcsPlusLoop;
            
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
        } else if (path.endsWith(".smi")){
            return new IteratingSMILESReader(
                    in, NoNotificationChemObjectBuilder.getInstance());
        } else {
            throw new Exception("Unrecognised filetype " + path);
        }

    }
    
}