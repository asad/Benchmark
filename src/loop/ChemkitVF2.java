package loop;

import org.openscience.cdk.interfaces.IMolecule;
import automorphism.vf2.AtomMapping;
import automorphism.vf2.VF2Automorphism;

public class ChemkitVF2 extends AbstractSubgraphIsomorphismLoop
        implements TimedSubgraphIsomorphismLoop {

    @Override
    public String getName() {
        return "ChemKitVF2";
    }

    @Override
    public void run(IMolecule query, IMolecule target) {
        VF2Automorphism matcher = new VF2Automorphism();
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
