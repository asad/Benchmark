package loop;

import gui.ImageGenerator;

import java.io.File;
import java.text.NumberFormat;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;

public abstract class AbstractSubgraphIsomorphismLoop {
    
    protected int numberOfResults;
    
    private long startTime;
    
    private long elapsedTime;
    
    public void startTimer() {
        startTime = System.currentTimeMillis();
    }
    
    public void stopTimer() {
        elapsedTime = System.currentTimeMillis() - startTime;
    }
    
    public long getTimeTaken() {
        return elapsedTime;
    }

    public int getResultCount() {
        return numberOfResults;
    }
    
    public void run(IMolecule query, List<IMolecule> targets) {
        startTimer();
        for (IMolecule target : targets) {
            run(query, target);
        }
        stopTimer();
    }
    
    public void generateImage(String outPutFileName, IAtomContainer query, IAtomContainer target, List<Map<Integer, Integer>> smsd) throws Exception {

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

    public abstract void run(IMolecule query, IMolecule target);
    
    public abstract String getName();
    
    public String toString() {
        return String.format("\t%st\t%d\t%ss\t%d",
                getName(), getTimeTaken(), getName(), getResultCount());
    }

}
