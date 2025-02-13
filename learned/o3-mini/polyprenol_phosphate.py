"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: polyprenol phosphate 
Definition: A prenol phosphate resulting from the formal condensation of the
terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
The algorithm verifies:
  1. The presence of a phosphate group attached (via an oxygen single bond) to
     an sp³-hybridized carbon that is allylic (i.e. has at least one double bond).
  2. The presence of a repeating isoprene-like motif ("C/C=C/CC"), indicative of a polyprenol chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    
    A polyprenol phosphate is defined as a compound derived from a polyprenol 
    (built of isoprene units) whose terminal allylic hydroxyl group has been esterified 
    with phosphoric acid. Here we search for an oxygen that connects a phosphorus (P) 
    (via a single bond) with an sp³ carbon that is allylic (i.e. adjacent to a C=C double bond).
    Additionally, we verify that an isoprene unit pattern ("C/C=C/CC") is present.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a polyprenol phosphate, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    found_phosphate_allylic = False
    # Loop over oxygen atoms (which generally form the bridge)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            neighbors = atom.GetNeighbors()
            # We must have at least one phosphorus and one carbon neighbor connected via single bonds.
            has_phosphorus = any(n.GetAtomicNum() == 15 for n in neighbors)
            has_carbon = any(n.GetAtomicNum() == 6 for n in neighbors)
            if has_phosphorus and has_carbon:
                # Verify that each connecting bond is a single bond.
                valid_bonds = True
                p_atom = None
                c_atom = None
                for nbr in neighbors:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                        valid_bonds = False
                        break
                    # Record one phosphorus and one carbon neighbor.
                    if nbr.GetAtomicNum() == 15 and p_atom is None:
                        p_atom = nbr
                    if nbr.GetAtomicNum() == 6 and c_atom is None:
                        c_atom = nbr
                if not valid_bonds or p_atom is None or c_atom is None:
                    continue
                # Check that the carbon is sp3 – expected for a former hydroxyl group.
                if c_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                # Check if this carbon is allylic: i.e., attached to at least one double bond with another carbon.
                allylic = False
                for bond_c in c_atom.GetBonds():
                    # Skip bonds connecting it back to the bridging oxygen.
                    other = bond_c.GetOtherAtom(c_atom)
                    if other.GetAtomicNum() == 6 and bond_c.GetBondType() == Chem.BondType.DOUBLE:
                        allylic = True
                        break
                if allylic:
                    found_phosphate_allylic = True
                    break  # We found a valid phosphate-allylic connection
    
    if not found_phosphate_allylic:
        return False, "No phosphate group found attached to an allylic carbon"
    
    # Additionally, verify the presence of an isoprene unit pattern.
    isoprene_pattern = Chem.MolFromSmarts("C/C=C/CC")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene unit pattern (C/C=C/CC) found indicative of a polyprenol chain"
    
    return True, "Contains phosphate group (via oxygen) attached to an allylic carbon and shows polyprenol repeating pattern"

# Example usage:
if __name__ == '__main__':
    # Replace test_smiles with any candidate SMILES for verification.
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\COP(O)(=O)OP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print(result, reason)