"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: polyprenol phosphate 
Definition: A prenol phosphate resulting from the formal condensation of the
terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

This algorithm uses two main criteria:
  1. The molecule contains at least one phosphate (a phosphorus atom) that is
     connected via an oxygen (by a single bond) to a carbon atom. This carbon must be
     sp³-hybridized and “allylic” (i.e. besides the O-bond it has at least one double bond
     to another carbon).
  2. The molecule displays an isoprene unit pattern (SMARTS: "C/C=C/CC") that is common
     to polyprenol chains.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    
    A polyprenol phosphate is defined as a compound derived from a polyprenol (built of isoprene
    units) whose terminal allylic hydroxyl group has been esterified with phosphoric acid.
    
    The algorithm checks that:
      1. There is at least one instance where a phosphorus atom is attached (via oxygen, single bond)
         to an sp3-hybridized carbon. For that carbon (which would have been the terminal allylic hydroxyl),
         at least one neighbor (other than the oxygen) is double-bonded to another carbon.
      2. The molecule contains an isoprene repeating pattern (SMARTS "C/C=C/CC").
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for a polyprenol phosphate.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    found_valid_phosphate = False
    # Instead of iterating over every oxygen, we iterate over all phosphorus atoms.
    for p in mol.GetAtoms():
        if p.GetAtomicNum() != 15:
            continue  # Not phosphorus
        # Loop over oxygens attached to phosphorus (via single bonds)
        for o in p.GetNeighbors():
            if o.GetAtomicNum() != 8:
                continue
            bp = mol.GetBondBetweenAtoms(p.GetIdx(), o.GetIdx())
            if bp is None or bp.GetBondType() != Chem.BondType.SINGLE:
                continue
            # For this oxygen, find a neighboring carbon (other than the phosphorus)
            for c in o.GetNeighbors():
                if c.GetAtomicNum() != 6:
                    continue
                # Ensure we are not going back to the phosphorus
                if c.GetIdx() == p.GetIdx():
                    continue
                bc = mol.GetBondBetweenAtoms(o.GetIdx(), c.GetIdx())
                if bc is None or bc.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Check that the carbon is sp3-hybridized (indicative of a former hydroxyl)
                if c.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                # Check the “allylic” nature of the carbon:
                # At least one other bond from c should be a double bond to another carbon.
                allylic = False
                for b in c.GetBonds():
                    # Skip the bond going back to oxygen
                    nb = b.GetOtherAtom(c)
                    if nb.GetIdx() == o.GetIdx():
                        continue
                    if nb.GetAtomicNum() != 6:
                        continue
                    if b.GetBondType() == Chem.BondType.DOUBLE:
                        allylic = True
                        break
                if allylic:
                    found_valid_phosphate = True
                    break  # Found a valid O-C(P) connection
            if found_valid_phosphate:
                break  # No need to check further oxygens for this phosphorus atom
        if found_valid_phosphate:
            break  # Found at least one valid phosphate connection
    
    if not found_valid_phosphate:
        return False, ("MISSED: No phosphate group found attached (via oxygen) to an sp3 carbon "
                       "that is allylic")
    
    # Check for isoprene repeating unit.
    # The SMARTS "C/C=C/CC" is used to identify a typical isoprene fragment.
    isoprene_pattern = Chem.MolFromSmarts("C/C=C/CC")
    if isoprene_pattern is None:
        return False, "Error in isoprene SMARTS pattern"
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, ("No isoprene unit pattern (C/C=C/CC) found indicative of a polyprenol chain")
    
    return True, ("Contains a phosphate (via oxygen) attached to an allylic sp3 carbon and exhibits "
                  "a polyprenol repeating pattern (isoprene units)")

# Example usage:
if __name__ == '__main__':
    # Replace 'test_smiles' below with any candidate SMILES string to test.
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\COP(O)(=O)OP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print(result, reason)