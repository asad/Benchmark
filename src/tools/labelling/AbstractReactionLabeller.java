package tools.labelling;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.Mapping;
import org.openscience.cdk.MoleculeSet;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.tools.manipulator.ReactionManipulator;

public class AbstractReactionLabeller {
    
    /**
     * A nasty hack necessary to get around a bug in the CDK
     */
    private boolean fixAtomMappingCastType = false;
    
    private void fixAtomMapping(IAtomContainer canonicalForm) {
        for (IAtom a : canonicalForm.atoms()) { 
            String v = (String) a.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
            if (v != null) {
                a.setProperty(
                        CDKConstants.ATOM_ATOM_MAPPING, Integer.valueOf(v));
            }
        }
    }
    
    protected void atomAtomMap(IMoleculeSet original, 
            IMoleculeSet clone, Map<IMolecule, int[]> permutationMap,
            Map<IAtom, IAtom> atomAtom) {
        for (int i = 0; i < original.getMoleculeCount(); ++i) {
            IMolecule mol = original.getMolecule(i);
            IMolecule mol2 = clone.getMolecule(i);
            int[] permutation = permutationMap.get(mol2);
            for (int j = 0; j < mol.getAtomCount(); ++j) {
                atomAtom.put(mol.getAtom(j), mol2.getAtom(permutation[j]));
            }
        }
    }
    
    protected Map<IAtom, IAtom> atomAtomMap(
            IReaction reaction, IReaction clone, 
            Map<IMolecule, int[]> permutationMap) {
        // create a Map of corresponding atoms for molecules 
        // (key: original Atom, value: clone Atom)
        Map<IAtom, IAtom> atomAtom = new Hashtable<IAtom, IAtom>();
        IMoleculeSet reactants = reaction.getReactants();
        IMoleculeSet clonedReactants = clone.getReactants();
        atomAtomMap(reactants, clonedReactants, permutationMap, atomAtom);
        
        IMoleculeSet products = reaction.getProducts();
        IMoleculeSet clonedProducts = clone.getProducts();
        atomAtomMap(products, clonedProducts, permutationMap, atomAtom);
        
        return atomAtom;
    }
    
    private void printAtomAtom(
            Map<IAtom, IAtom> atomAtom, IReaction reaction, IReaction clone) {
        for (IAtom key : atomAtom.keySet()) {
            IAtomContainer keyAC = 
                ReactionManipulator.getRelevantAtomContainer(reaction, key);
            int keyIndex = keyAC.getAtomNumber(key);
            IAtom value = atomAtom.get(key);
            IAtomContainer valueAC = 
                ReactionManipulator.getRelevantAtomContainer(clone, value);
            int valueIndex = valueAC.getAtomNumber(value);
            System.out.println(
                    "key " + keyIndex + key.getSymbol() 
                    + " mapped to " + valueIndex + value.getSymbol());
        }
    }
    
    protected List<IMapping> cloneMappings(
            IReaction reaction, Map<IAtom, IAtom> atomAtomMap) {
        // clone the mappings
        int numberOfMappings = reaction.getMappingCount();
        List<IMapping> map = new ArrayList<IMapping>();
        for (int mappingIndex = 0; mappingIndex < numberOfMappings; mappingIndex++) {
            IMapping mapping = reaction.getMapping(mappingIndex);
            map.add(cloneMapping(mapping, atomAtomMap));
        }
        return map;
    }
    
    protected IMapping cloneMapping(
            IMapping mapping, Map<IAtom, IAtom> atomAtomMap) {
        IChemObject keyChemObj0 = mapping.getChemObject(0);
        IChemObject keyChemObj1 = mapping.getChemObject(1);
        IChemObject co0 = (IChemObject) atomAtomMap.get(keyChemObj0);
        IChemObject co1 = (IChemObject) atomAtomMap.get(keyChemObj1);
        
        // THIS IS STUPID : BLAME THE IDIOT WHO FAILED TO PUT SET METHODS IN
        // IMAPPING (OR IREACTION, FOR THAT MATTER)
        if (co0 == null) {
            co0 = new ChemObject();
        }
        if (co1 == null) {
            co1 = new ChemObject();
        }
        return new Mapping(co0, co1);
    }
    
    protected Map<IChemObject, Integer> makeIndexMap(IReaction reaction) {
        Map<IChemObject, Integer> indexMap = 
            new HashMap<IChemObject, Integer>();
        List<IAtomContainer> all = 
            ReactionManipulator.getAllAtomContainers(reaction);
        int globalIndex = 0;
        for (IAtomContainer ac : all) { 
            for (IAtom atom : ac.atoms()) {
                indexMap.put(atom, globalIndex);
                globalIndex++;
            }
        }
        return indexMap;
    }
    
    /**
     * Clone and Sort the mappings based on the order of the first object 
     * in the mapping (which is assumed to be the reactant).
     * 
     * @param reaction
     */
    protected void cloneAndSortMappings(
            IReaction reaction, IReaction copyOfReaction, 
            Map<IMolecule, int[]> permutationMap) {
        
        Map<IAtom, IAtom> atomAtomMap = atomAtomMap(
                reaction, copyOfReaction, permutationMap);
        printAtomAtom(atomAtomMap, reaction, copyOfReaction);
        List<IMapping> map = cloneMappings(reaction, atomAtomMap);
        sortMappings(copyOfReaction, map);
    }
    
    protected void sortMappings(IReaction reaction, List<IMapping> map) {
        // make a lookup for the indices of the atoms 
        final Map<IChemObject, Integer> indexMap = makeIndexMap(reaction);
        Comparator<IMapping> mappingSorter = new Comparator<IMapping>() {

            @Override
            public int compare(IMapping o1, IMapping o2) {
                IChemObject o10 = o1.getChemObject(0);
                IChemObject o20 = o2.getChemObject(0);
                if (o20 == null || o10 == null) return 0;
                Integer o10i = indexMap.get(o10); 
                Integer o20i = indexMap.get(o20);
                if (o10i == null || o20i == null) return 0;
                return o10i.compareTo(o20i);
            }
            
        };
        Collections.sort(map, mappingSorter);
        int mappingIndex = 0;
        for (IMapping mapping : map) {
            IChemObject o0 = mapping.getChemObject(0);
            if (o0 != null) {
                System.out.println("setting " + mappingIndex);
                o0.setProperty(CDKConstants.ATOM_ATOM_MAPPING, mappingIndex);
            } 
            IChemObject o1 = mapping.getChemObject(1);
            if (o1 != null) {
                o1.setProperty(CDKConstants.ATOM_ATOM_MAPPING, mappingIndex);
                System.out.println("setting " + mappingIndex);
            }
            reaction.addMapping(mapping);
            printMapping(mapping, reaction);
            mappingIndex++;
        }
    }
    
    private void printMapping(IMapping mapping, IReaction reaction) {
        IChemObject c0 = (IChemObject) mapping.getChemObject(0);
        int a0i = -1;
        if (c0 != null && c0 instanceof IAtom) {
            IAtom a0 = (IAtom) c0;
            IAtomContainer ac = 
                ReactionManipulator.getRelevantAtomContainer(reaction, a0);
            if (ac != null) {
                a0i = ac.getAtomNumber(a0);
            }
        }
        IChemObject c1 = (IChemObject) mapping.getChemObject(1);
        int a1i = -1;
        if (c1 != null && c1 instanceof IAtom) {
            IAtom a1 = (IAtom) c1;
            IAtomContainer ac = 
                ReactionManipulator.getRelevantAtomContainer(reaction, a1);
            if (ac != null) {
                a1i = ac.getAtomNumber(a1);
            }
        }
        System.out.println("(" + a0i + ", " + a1i + ")");
    }
    
    protected IMoleculeSet canoniseMoleculeSet(IMoleculeSet moleculeSet,
            ICanonicalMoleculeLabeller labeller,
            Map<IMolecule, int[]> permutationMap) {
        IMoleculeSet canonicalMolecules = new MoleculeSet();
        for (IAtomContainer atomContainer : moleculeSet.atomContainers()) {
            IAtomContainer canonicalForm = 
                labeller.getCanonicalMolecule(atomContainer);
            if (fixAtomMappingCastType) { fixAtomMapping(canonicalForm); }
            IMolecule canonicalMolecule = 
                canonicalForm.getBuilder().newInstance(
                        IMolecule.class, canonicalForm); 
            permutationMap.put(
                    canonicalMolecule, 
                    labeller.getCanonicalPermutation(atomContainer));
            canonicalMolecules.addMolecule(canonicalMolecule);
        }
        return canonicalMolecules;
    }
    
    public IReaction labelReaction(
            IReaction reaction, ICanonicalMoleculeLabeller labeller) {
        System.out.println("labelling");
        IReaction canonReaction = new Reaction();
        
        Map<IMolecule, int[]> permutationMap = new HashMap<IMolecule, int[]>();
        
        canonReaction.setProducts(
                canoniseMoleculeSet(
                        reaction.getProducts(), labeller, permutationMap));
        canonReaction.setReactants(
                canoniseMoleculeSet(
                        reaction.getReactants(), labeller, permutationMap));
        
        cloneAndSortMappings(reaction, canonReaction, permutationMap);
        return canonReaction;
    }

}
