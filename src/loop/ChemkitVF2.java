package loop;

import chemkit.substructure.VF2;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;

public class ChemkitVF2 extends AbstractSubgraphIsomorphismLoop
        implements TimedSubgraphIsomorphismLoop {

    @Override
    public String getName() {
        return "ChemKitVF2";
    }

    /**
     * 
     * @param query
     * @param target
     */
    @Override
    public void run(IMolecule query, IMolecule target) {
        if (query.getAtomCount() <= target.getAtomCount()) {
            VF2 matcher = new VF2(true, false);
            try {
                matcher.set(query, target);
            } catch (CDKException ex) {
                Logger.getLogger(ChemkitVF2.class.getName()).log(Level.SEVERE, null, ex);
            }
            if (matcher.isSubgraph()) {
                numberOfResults++;
            }
        }
    }
}
