package helper;

import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.labelling.CanonicalLabellingAdaptor;
import org.openscience.smsd.labelling.ICanonicalMoleculeLabeller;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

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
     */
    public String getCanonicalSMILES(IAtomContainer mol) {
        this.molecule = ExtAtomContainerManipulator.makeDeepCopy(mol);
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