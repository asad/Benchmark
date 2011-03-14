package isomorphism.vf2.bond;

import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * 
 * @author Asad
 */
public class BondMapping implements Cloneable {

    private IAtomContainer a;
    private IAtomContainer b;
    private Map<IBond, IBond> mapping;

    /**
     * 
     * @param a source mol
     * @param b target mol
     */
    public BondMapping(IAtomContainer a, IAtomContainer b) {
        this.a = a;
        this.b = b;
        this.mapping = new HashMap<IBond, IBond>();
    }

    /**
     * 
     * @param atom1
     * @param atom2
     */
    public void add(IBond atom1, IBond atom2) {
        mapping.put(atom1, atom2);
    }

    @Override
    public String toString() {
        String s = "[";
        for (IBond key : mapping.keySet()) {
            int keyIndex = a.getBondNumber(key);
            int valueIndex = b.getBondNumber(mapping.get(key));
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
     *size of the mapping
     * @return size of the mapping
     */
    public int size() {
        return mapping.size();
    }

    boolean containsQueryBond(IBond atom) {
        return mapping.containsKey(atom);
    }

    Iterable<IBond> queryBonds() {
        return mapping.keySet();
    }

    IBond getMappedTargetBond(IBond atom) {
        return mapping.get(atom);
    }

    boolean containsTargetBond(IBond atom) {
        return mapping.containsValue(atom);
    }

    Map<IBond, IBond> getMapping() {
        return mapping;
    }
}