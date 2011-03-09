package loop;

import java.util.List;

import org.openscience.cdk.interfaces.IMolecule;

public interface TimedSubgraphIsomorphismLoop {

    public void run(IMolecule query, List<IMolecule> targets);
    
    public int getResultCount();
}
