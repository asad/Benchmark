package loop;

import java.util.List;
import java.util.Map;

import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.smsd.algorithm.mcsplus.MCSPlusHandler;
import org.openscience.smsd.tools.MolHandler;

public class MCSPlusLoop extends AbstractSubgraphIsomorphismLoop implements
        TimedSubgraphIsomorphismLoop {

    public String getName() {
        return "MCSPlus";
    }
    
    @Override
    public void run(IMolecule query, IMolecule target) {
        MCSPlusHandler mcsplus = new MCSPlusHandler();
        mcsplus.set(new MolHandler(query, false, false), new MolHandler(target, false, false));
        mcsplus.searchMCS(true);

        if (mcsplus.getFirstMapping() != null && !mcsplus.getFirstMapping().isEmpty()) {
//            List<Map<Integer, Integer>> map = mcsplus.getAllMapping();
            numberOfResults++;
        }
    }

}
