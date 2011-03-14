package smsd.vf2.atom;

import java.util.HashMap;
import java.util.Map;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * 
 * @author Asad
 */
public class AtomMapping implements Cloneable {

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
     *size of the mapping
     * @return size of the mapping
     */
    public int size() {
        return mapping.size();
    }

    boolean containsQueryAtom(IAtom atom) {
        return mapping.containsKey(atom);
    }

    Iterable<IAtom> queryAtoms() {
        return mapping.keySet();
    }

    IAtom getMappedTargetAtom(IAtom atom) {
        return mapping.get(atom);
    }

    boolean containsTargetAtom(IAtom atom) {
        return mapping.containsValue(atom);
    }

    Map<IAtom, IAtom> getMapping() {
        return mapping;
    }
}