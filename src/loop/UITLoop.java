package loop;

import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;

public class UITLoop extends AbstractSubgraphIsomorphismLoop
                     implements TimedSubgraphIsomorphismLoop {
    
    public String getName() {
        return "UIT";
    }
   
    public void run(IMolecule query, IMolecule target) {
        try {
            List bondMapping = UniversalIsomorphismTester.getSubgraphMaps(target, query);
            List<List<RMap>> sol = UniversalIsomorphismTester.makeAtomsMapsOfBondsMaps(bondMapping, target, query);
            if (sol.size() > 0) {
                numberOfResults++;
            }
        } catch (CDKException cdke) {
            
        }
    }

    @Override
    public int getResultCount() {
        return numberOfResults;
    }

}
