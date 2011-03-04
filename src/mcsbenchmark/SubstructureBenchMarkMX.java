package mcsbenchmark;

import com.metamolecular.mx.io.mdl.SDFileReader;
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

/**
 *
 * @author Asad
 */
public class SubstructureBenchMarkMX {

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

        SDFileReader qFileReader = read(queryFilePath);

        if (qFileReader == null) {
            throw new IOException("Unknown input type ");
        }
        List<Molecule> targets = new ArrayList<Molecule>();
        SDFileReader tFileReader = read(targetFilePath);
        while (tFileReader.hasNextRecord()) {
            tFileReader.nextRecord();
            targets.add(tFileReader.getMolecule());
        }

        int counter = 0;
        while (qFileReader.hasNextRecord()) {
            qFileReader.nextRecord();
            Molecule query = qFileReader.getMolecule();
            int mxSolutionCount = 0;
            long tMX0 = System.currentTimeMillis();
            for (Molecule target : targets) {
                mxSolutionCount += getMXSolutionCount(query, target);
            }
            long timeNow = System.currentTimeMillis();
            long mxTime = timeNow - tMX0;

            String out = String.format("%d MXt %d MXs %d",
                    counter, mxTime, mxSolutionCount);
            System.out.println(out);
            counter++;
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

    /**
     *
     * @param file
     * @return
     * @throws FileNotFoundException
     * @throws IOException  
     */
    public static SDFileReader read(String file) throws FileNotFoundException, IOException {

        FileReader in = new FileReader(file);
        return new SDFileReader(file);

    }
}
