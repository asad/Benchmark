package loop;

import org.openscience.cdk.interfaces.IMolecule;
import smsd.algorithm.mcsplus.MCSPlusHandler;

public class MCSPlusLoop extends AbstractSubgraphIsomorphismLoop implements
        TimedSubgraphIsomorphismLoop {

    @Override
    public String getName() {
        return "MCSPlus";
    }

    /**
     * 
     * @param query
     * @param target
     */
    @Override
    public void run(IMolecule query, IMolecule target) {
        if (query.getAtomCount() <= target.getAtomCount()) {
            MCSPlusHandler mcsplus = new MCSPlusHandler();
            mcsplus.set(query, target);
            mcsplus.searchMCS(true, true);

            if (mcsplus.getFirstMapping() != null && !mcsplus.getFirstMapping().isEmpty()) {
                numberOfResults++;
            }
        }
    }
}
