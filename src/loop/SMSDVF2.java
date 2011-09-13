package loop;

import org.openscience.cdk.interfaces.IMolecule;
import smsd.algorithm.vflib.VF2Sub;

public class SMSDVF2 extends AbstractSubgraphIsomorphismLoop
        implements TimedSubgraphIsomorphismLoop {

    @Override
    public String getName() {
        return "SMSD";
    }

    /**
     * 
     * @param query
     * @param target
     */
    @Override
    public void run(IMolecule query, IMolecule target) {
        if (query.getAtomCount() <= target.getAtomCount()) {
            VF2Sub matcher = new VF2Sub(true, true);
            matcher.set(query, target);
            if (!matcher.isSubgraph()) {
                numberOfResults++;
            }
        }
    }
}
