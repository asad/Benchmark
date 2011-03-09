package loop;

import mcsbenchmark.VF2;
import mcsbenchmark.VF2.AtomMapping;

import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smiles.SmilesGenerator;

public class VFLibLoop extends AbstractSubgraphIsomorphismLoop 
                       implements TimedSubgraphIsomorphismLoop {
    
    public String getName() {
        return "VFLib";
    }
    
    public void run(IMolecule query, IMolecule target) {
        SmilesGenerator smilesGenerator = new SmilesGenerator();
        VF2 matcher = new VF2();
        AtomMapping mapping = matcher.isomorphism(query, target);
//        List<AtomMapping> mappings = matcher.getAllMappings(query, target);
        if (!mapping.isEmpty()) {
//            System.out.println(mapping);
//            System.out.println(smilesGenerator.createSMILES(query));
//            System.out.println(smilesGenerator.createSMILES(target));
            numberOfResults++;
        }
//        return mappings.size();
    }
}
