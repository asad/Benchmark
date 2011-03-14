package loop;

import org.openscience.cdk.interfaces.IMolecule;
import chemkit.vf2.AtomMapping;
import chemkit.vf2.VF2;

public class ChemkitVF2 extends AbstractSubgraphIsomorphismLoop
        implements TimedSubgraphIsomorphismLoop {

    @Override
    public String getName() {
        return "ChemKitVF2";
    }

    @Override
    public void run(IMolecule query, IMolecule target) {
        VF2 matcher = new VF2();
        AtomMapping mapping = matcher.isomorphism(query, target);
        if (!mapping.isEmpty()) {
            numberOfResults++;
        }
//        VFMapper matcher = new VFMapper(query);
//        if (matcher.hasMap(target)) {
//            numberOfResults++;
//        }
    }
}
