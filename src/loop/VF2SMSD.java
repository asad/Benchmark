package loop;

import org.openscience.cdk.interfaces.IMolecule;
import vf2.old.VFMapper;

public class VF2SMSD extends AbstractSubgraphIsomorphismLoop
        implements TimedSubgraphIsomorphismLoop {

    @Override
    public String getName() {
        return "SMSD";
    }

    @Override
    public void run(IMolecule query, IMolecule target) {
        if (query.getAtomCount() <= target.getAtomCount()) {
            VFMapper matcher = new VFMapper(query);
            if (matcher.getFirstMap(target).isEmpty()) {
                numberOfResults++;
            }
        }
    }
}
