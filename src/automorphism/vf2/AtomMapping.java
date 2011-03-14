package automorphism.vf2;

import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * 
 * @author Asad
 */
public class AtomMapping {

    private IAtomContainer a;
    private IAtomContainer b;
    private Map<IAtom, IAtom> mapping;

    /**
     * 
     * @param a source mol
     * @param b target mol
     */
    public AtomMapping(IAtomContainer a, IAtomContainer b) {
        this.a = a;
        this.b = b;
        this.mapping = new HashMap<IAtom, IAtom>();
    }

    /**
     * 
     * @param atom1
     * @param atom2
     */
    public void add(IAtom atom1, IAtom atom2) {
        mapping.put(atom1, atom2);
    }

    @Override
    public String toString() {
        String s = "[";
        for (IAtom key : mapping.keySet()) {
            int keyIndex = a.getAtomNumber(key);
            int valueIndex = b.getAtomNumber(mapping.get(key));
            s += keyIndex + ":" + valueIndex + "|";
        }
        return s + "]";
    }

    /**
     * 
     * @return true if 'a' is not a subgraph of 'b'
     */
    public boolean isEmpty() {
        return mapping.isEmpty();
    }

    /**
     * 
     * clear mapping
     */
    public void clear() {
        mapping.clear();
    }
    
    /**
     * 
     * mapping size
     * @return 
     */
    public int mappingCount() {
        return mapping.size();
    }
}