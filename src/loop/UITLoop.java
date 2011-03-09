package loop;

import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.smiles.SmilesGenerator;

public class UITLoop extends AbstractSubgraphIsomorphismLoop
                     implements TimedSubgraphIsomorphismLoop {
    
    public String getName() {
        return "UIT";
    }
   
    public void run(IMolecule query, IMolecule target) {
        try {
            List bondMapping = UniversalIsomorphismTester.getSubgraphMaps(target, query);
            List<List<RMap>> sol = UniversalIsomorphismTester.makeAtomsMapsOfBondsMaps(bondMapping, target, query);
            //         List<IAtomContainer> list = UniversalIsomorphismTester.getOverlaps(query, target);
    
            //         CDKMCSHandler cdkmcs = new CDKMCSHandler();
            //         cdkmcs.set(new MolHandler(query, false, false), new MolHandler(target, false, false));
            //         cdkmcs.searchMCS(true);
            //
            //         if (cdkmcs.getFirstMapping() != null && !cdkmcs.getFirstMapping().isEmpty()) {
            //             List<Map<Integer, Integer>> map = cdkmcs.getAllMapping();
            //             if (mapSol.containsKey(counter)) {
            //                 int solcount = map.size();
            //                 int solsize = cdkmcs.getFirstAtomMapping().size();
            //                 List<Integer> l = mapSol.get(counter);
            //                 l.set(2, solcount);
            //                 l.set(3, solsize);
            //                 mapSol.put(counter, l);
            //             }
            ////             generateImage(counter + "_UIT", query, target, map);
            //             return 1;
            //         } else {
            //             return 0;
            //         }
            if (sol.size() > 0) {
//                SmilesGenerator smilesGenerator = new SmilesGenerator();
                //            System.out.println(sol);
//                System.out.println(smilesGenerator.createSMILES(query));
//                System.out.println(smilesGenerator.createSMILES(target));
    
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
