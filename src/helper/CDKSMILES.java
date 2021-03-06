package helper;

import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import tools.labelling.CanonicalLabellingAdaptor;
import tools.labelling.ICanonicalMoleculeLabeller;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK. e-mail: asad@ebi.ac.uk
 */
public class CDKSMILES {

    private final ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(CDKSMILES.class);
    private IAtomContainer molecule;

    public CDKSMILES() {
    }

    /**
     * 
     * @param mol
     * @return
     * @throws CloneNotSupportedException  
     */
    public String getCanonicalSMILES(IAtomContainer mol) throws CloneNotSupportedException {
        this.molecule = (GraphAtomContainer) mol.clone();
        if (!isPseudoAtoms()) {
            ICanonicalMoleculeLabeller canonLabeler = new CanonicalLabellingAdaptor();
            molecule = canonLabeler.getCanonicalMolecule(molecule);
        }
        SmilesGenerator sg = new SmilesGenerator(true);
        String smiles = "";
        try {
            smiles = sg.createSMILESWithoutCheckForMultipleMolecules(molecule, false, new boolean[molecule.getBondCount()]);
        } catch (CDKException ex) {
            logger.warn("SmilesGenerator error: ", ex.getMessage());
        }
        return smiles;
    }

    /**
     * 
     * @param mol
     * @return
     */
    public String getSMILES(IAtomContainer mol) {
        SmilesGenerator sg = new SmilesGenerator(true);
        String smiles = "";
        try {
            smiles = sg.createSMILESWithoutCheckForMultipleMolecules(molecule, false, new boolean[molecule.getBondCount()]);
        } catch (CDKException ex) {
            logger.warn("SmilesGenerator error: ", ex.getMessage());
        }
        return smiles;
    }

    private boolean isPseudoAtoms() {
        for (IAtom atoms : molecule.atoms()) {
            if (atoms instanceof IPseudoAtom || atoms instanceof PseudoAtom) {
                return true;
            }
        }
        return false;
    }
}