package loop;

import isomorphism.vf2.atom.VFAtomMapper;
import org.openscience.cdk.interfaces.IMolecule;

public class SMSDVF2 extends AbstractSubgraphIsomorphismLoop
        implements TimedSubgraphIsomorphismLoop {

    @Override
    public String getName() {
        return "SMSD";
    }

    @Override
    public void run(IMolecule query, IMolecule target) {
        if (query.getAtomCount() <= target.getAtomCount()) {
            VFAtomMapper matcher = new VFAtomMapper(query);
            if (!matcher.getFirstMap(target).isEmpty()) {
                numberOfResults++;
            }
//            VFBondMapper matcher = new VFBondMapper(query);
//            if (!matcher.getFirstMap(target).isEmpty()) {
//                numberOfResults++;
//            }
        }
    }
}
