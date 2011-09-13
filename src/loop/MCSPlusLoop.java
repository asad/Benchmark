package loop;

import org.openscience.cdk.interfaces.IMolecule;
import smsd.algorithm.mcsplus.MCSPlusHandler;

public class MCSPlusLoop extends AbstractSubgraphIsomorphismLoop implements
        TimedSubgraphIsomorphismLoop {

    @Override
    public String getName() {
        return "MCSPlus";
    }

    @Override
    public void run(IMolecule query, IMolecule target) {
        MCSPlusHandler mcsplus = new MCSPlusHandler();
        mcsplus.set(query, target);
        mcsplus.searchMCS(true, true);

        if (mcsplus.getFirstMapping() != null && !mcsplus.getFirstMapping().isEmpty()) {
            numberOfResults++;
        }
    }
}
