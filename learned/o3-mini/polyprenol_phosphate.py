"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: polyprenol phosphate 
Definition: A prenol phosphate resulting from the formal condensation of the
terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
The algorithm verifies:
  1. The presence of a phosphate (P) atom whose oxygen (via a single bond)
     is attached to an sp³ carbon that is allylic (i.e. has a neighboring double bond).
  2. A repeating isoprene-like motif, represented by the SMARTS "C/C=C/CC",
     which is expected for polyprenol chains.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    
    A polyprenol phosphate is defined as a compound derived from a polyprenol
    (a series of isoprene units) whose terminal allylic hydroxyl group has been esterified
    with phosphoric acid. This function first looks for a phosphate group connected via
    an oxygen to a carbon that is allylic (i.e. adjacent to a C=C double bond). Second,
    it searches for an isoprene-like pattern ("C/C=C/CC") as evidence of a polyprenol chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a polyprenol phosphate, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flag to mark whether we have found a phosphate group attached to an allylic carbon.
    found_phosphate_allylic = False
    
    # Loop through atoms to locate a phosphorus atom.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Phosphorus
            # Inspect neighboring atoms of phosphorus.
            for nbr in atom.GetNeighbors():
                # We are interested in oxygen neighbors.
                if nbr.GetAtomicNum() == 8:
                    # Get the bond between phosphorus and oxygen.
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    # Looking for a single bond (i.e. not a P=O bond).
                    if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                        # Now, check that this oxygen links to an organic (carbon) group.
                        for o_nbr in nbr.GetNeighbors():
                            if o_nbr.GetIdx() == atom.GetIdx():
                                continue  # Skip the phosphorus already considered.
                            if o_nbr.GetAtomicNum() == 6:  # Carbon neighbor
                                # Check that the candidate carbon is sp3,
                                # as expected for a former hydroxyl-bearing carbon.
                                if o_nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                                    # Now determine if this carbon is allylic:
                                    # i.e. at least one bond from this carbon is a double bond (C=C).
                                    allylic = False
                                    for bond_c in o_nbr.GetBonds():
                                        # Identify the atom at the other end of the bond.
                                        other = bond_c.GetOtherAtom(o_nbr)
                                        if other.GetAtomicNum() == 6 and bond_c.GetBondType() == Chem.BondType.DOUBLE:
                                            allylic = True
                                            break
                                    if allylic:
                                        found_phosphate_allylic = True
                                        break
                        if found_phosphate_allylic:
                            break
            if found_phosphate_allylic:
                break
                
    if not found_phosphate_allylic:
        return False, "No phosphate group found attached to an allylic carbon"
    
    # Additionally, search for an isoprene unit pattern as a proxy for a polyprenol chain.
    # The SMARTS "C/C=C/CC" represents a basic isoprene fragment.
    isoprene_pattern = Chem.MolFromSmarts("C/C=C/CC")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene unit pattern (C/C=C/CC) found indicative of a polyprenol chain"
    
    return True, "Contains phosphate group esterified to an allylic carbon and shows isoprene repeating pattern indicative of a polyprenol chain"

# Example usage:
if __name__ == '__main__':
    # A sample SMILES for testing. Replace with any of the provided examples.
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\COP(O)(=O)OP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print(result, reason)