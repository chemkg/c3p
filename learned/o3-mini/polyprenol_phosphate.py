"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: polyprenol phosphate 
Definition: A prenol phosphate resulting from the formal condensation of the
terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
The algorithm verifies:
  1. That there is at least one oxygen atom that forms a single bond with a phosphorus atom
     and also connects to an sp³-hybridized carbon that is allylic (i.e. the carbon has at least one
     double bond connected to another carbon).
  2. The molecule contains an isoprene unit pattern (e.g., "C/C=C/CC") indicative of a polyprenol chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    
    A polyprenol phosphate is defined as a compound derived from a polyprenol (built of isoprene units)
    whose terminal allylic hydroxyl group has been esterified with phosphoric acid (or a related diphosphate).
    
    The algorithm checks that:
      1. There is an oxygen that is connected by single bonds to at least one phosphorus atom and at least
         one carbon atom. For each candidate carbon attached via such an oxygen, we require that it is sp3-hybridized
         (as expected for a former hydroxyl) and that it is “allylic” (has at least one double bond to another carbon).
      2. The molecule contains an isoprene repeating unit pattern.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a polyprenol phosphate, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    found_valid_connection = False
    # Loop over every oxygen atom
    for oxygen in mol.GetAtoms():
        if oxygen.GetAtomicNum() != 8:
            continue
        # Look at each neighbor of this oxygen.
        # We want at least one phosphorus neighbor (via a single bond) and at least one sp3 carbon neighbor (via a single bond)
        phosphorus_found = any(n.GetAtomicNum() == 15 and mol.GetBondBetweenAtoms(oxygen.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.SINGLE 
                               for n in oxygen.GetNeighbors())
        if not phosphorus_found:
            continue  # no P attached, try next oxygen
        # Now check each carbon neighbor attached via a single bond.
        for n in oxygen.GetNeighbors():
            if n.GetAtomicNum() != 6:
                continue
            bond = mol.GetBondBetweenAtoms(oxygen.GetIdx(), n.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Check that the carbon is sp3-hybridized (as expected for a former hydroxyl carbon)
            if n.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # Check if this carbon is allylic: that is, if any of its bonds is a double bond to another carbon.
            allylic = False
            for b in n.GetBonds():
                # We skip the bond that goes back to the oxygen we already considered.
                if b.GetOtherAtom(n).GetIdx() == oxygen.GetIdx():
                    continue
                # Check for a double bond to another carbon.
                if b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(n).GetAtomicNum() == 6:
                    allylic = True
                    break
            if allylic:
                found_valid_connection = True
                break  # found at least one valid O-C(P) connection
        if found_valid_connection:
            break  # no need to check further oxygens
    
    if not found_valid_connection:
        return False, "No phosphate group found attached (via oxygen) to an sp3 carbon that is allylic"
    
    # Now verify the presence of an isoprene-like pattern.
    # The SMARTS "C/C=C/CC" should match a typical isoprene unit.
    isoprene_pattern = Chem.MolFromSmarts("C/C=C/CC")
    if isoprene_pattern is None:
        return False, "Error in isoprene SMARTS pattern"
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene unit pattern (C/C=C/CC) found indicative of a polyprenol chain"
    
    return True, "Contains phosphate (via oxygen) attached to an allylic sp3 carbon and exhibits polyprenol repeating pattern"

# Example usage:
if __name__ == '__main__':
    # Replace test_smiles below with any candidate SMILES string to test.
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\COP(O)(=O)OP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print(result, reason)