package loop;

import chemkit.AtomMapping;
import chemkit.VF2;
import org.openscience.cdk.interfaces.IMolecule;

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
    }
}
