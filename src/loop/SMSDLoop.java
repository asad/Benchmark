package loop;

import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.smsd.algorithm.vflib.VFlibMCSHandler;
import org.openscience.smsd.tools.MolHandler;

public class SMSDLoop extends AbstractSubgraphIsomorphismLoop
                      implements TimedSubgraphIsomorphismLoop {

    public String getName() {
        return "SMSD";
    }
    
    public void run(IMolecule query, IMolecule target) {
        VFlibMCSHandler mcs = null;
        mcs = new VFlibMCSHandler();
        mcs.set(new MolHandler(query, false, false), new MolHandler(target, false, false));
        mcs.searchMCS(true);
        if (mcs.getFirstMapping() != null && !mcs.getFirstMapping().isEmpty()) {
            numberOfResults++;
        }
    }
}
