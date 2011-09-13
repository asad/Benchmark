package loop;

import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.SMILESWriter;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;

public class UITLoop extends AbstractSubgraphIsomorphismLoop
        implements TimedSubgraphIsomorphismLoop {

    @Override
    public String getName() {
        return "UIT";
    }

    @Override
    public void run(IMolecule query, IMolecule target) {
        SMILESWriter sw=new SMILESWriter();
        if (query.getAtomCount() <= target.getAtomCount()) {
            try {
                List bondMapping = UniversalIsomorphismTester.getSubgraphMap(target, query);
                List<RMap> sol = UniversalIsomorphismTester.makeAtomsMapOfBondsMap(bondMapping, target, query);
                if (sol != null && sol.size() > 0) {
                    sw.writeMolecule(query);
                    numberOfResults++;
                }
            } catch (CDKException cdke) {
            }
        }
    }

    @Override
    public int getResultCount() {
        return numberOfResults;
    }
}
