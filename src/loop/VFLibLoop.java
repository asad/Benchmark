package loop;

import mcsbenchmark.VF2;
import mcsbenchmark.VF2.AtomMapping;

import org.openscience.cdk.interfaces.IMolecule;

public class VFLibLoop extends AbstractSubgraphIsomorphismLoop 
                       implements TimedSubgraphIsomorphismLoop {
    
    public String getName() {
        return "VFLib";
    }
    
    public void run(IMolecule query, IMolecule target) {
        VF2 matcher = new VF2();
        AtomMapping mapping = matcher.isomorphism(query, target);
        if (!mapping.isEmpty()) {
            numberOfResults++;
        }
    }
}
