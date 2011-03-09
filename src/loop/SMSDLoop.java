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
//            List<Map<Integer, Integer>> map = mcs.getAllMapping();
//            if (mapSol.containsKey(counter)) {
//                int solcount = map.size();
//                int solsize = mcs.getFirstAtomMapping().size();
//                List<Integer> l = mapSol.get(counter);
//                l.set(0, solcount);
//                l.set(1, solsize);
//                mapSol.put(counter, l);
//            }
//            generateImage(counter + "_SMSD", query, target, map);
            numberOfResults++;
        }
    }
}
