"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
Definition: A carbohydrate derivative in which one or more of the oxygens or hydroxy groups 
of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen or any group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    
    The function uses a heuristic approach:
      1. The molecule must contain at least one sulfur atom.
      2. The molecule is searched for a candidate sugar ring (ring of size 5 or 6)
         that has a majority of carbons plus one heteroatom (oxygen or sulfur) and several hydroxyl groups.
      3. If the candidate sugar ring contains a ring heteroatom that is sulfur, this
         suggests the classical thiosugar (a thiopyranose or thiofuranose).
      4. Alternatively, if one of the sugar ring carbons is substituted by an S atom 
         (outside the ring) and the ring carries several -OH groups, we interpret this as a thiosugar derivative.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str): Tuple where bool is True if molecule is classified as a thiosugar,
                     and str provides the classification rationale.
        
    If the SMILES cannot be parsed or no suitable features are found, returns (False, reason).
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens for correct -OH group detection
    molH = Chem.AddHs(mol)
    
    # Requirement 1: must have at least one sulfur atom.
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atom present, so not a thiosugar"
        
    # Retrieve ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Loop over rings to find one that could be a sugar ring candidate.
    for ring in rings:
        # Consider only rings of typical sugar sizes (5 or 6 members)
        if len(ring) not in [5, 6]:
            continue
            
        # Count carbons in the ring and gather heteroatoms (O and S)
        num_carbons = 0
        hetero_atoms = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                num_carbons += 1
            elif atom.GetAtomicNum() in (8, 16):  # oxygen or sulfur
                hetero_atoms.append(atom)
                
        # A typical sugar ring usually has at least 4 carbons and one ring heteroatom.
        if num_carbons < 4 or len(hetero_atoms) < 1:
            continue

        # Case 1: the ring heteroatom (normally oxygen in sugars) is replaced by sulfur.
        if len(hetero_atoms) == 1 and hetero_atoms[0].GetAtomicNum() == 16:
            return True, "Thiosugar identified: sugar ring with sulfur replacing the ring oxygen."

        # Case 2: the ring appears to be a carbohydrate (many hydroxyl groups on ring substituents)
        # and one of the ring carbons is substituted by a sulfur (i.e. -S or -SR) replacing an -OH.
        # We count substituent –OH groups attached to ring atoms as typical for sugars.
        oh_count = 0
        for idx in ring:
            ring_atom = molH.GetAtomWithIdx(idx)
            # Look for neighbors outside the ring that are oxygen and bound to at least one hydrogen.
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms within the ring
                if nbr.GetSymbol() == "O":
                    # Check if this oxygen has at least one hydrogen neighbor (an -OH group).
                    if any(nei.GetSymbol() == "H" for nei in nbr.GetNeighbors()):
                        oh_count += 1
        
        # Now check if any ring carbon has an external sulfur substituent.
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            if ring_atom.GetAtomicNum() != 6:
                continue  # only consider carbons in the ring
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 16:
                    # We require that the ring has a fair number of hydroxyl groups to be sugar-like.
                    if oh_count >= 2:  # arbitrary threshold for a sugar-like pattern
                        return True, "Thiosugar identified: sugar ring with a sulfur substituent replacing a hydroxyl group."
                        
    return False, "No thiosugar substructure found"
    
# Example usage:
if __name__ == '__main__':
    # Test with one thiosugar example: 6-thio-beta-D-galactose
    test_smiles = "[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CS)O)O)O)O"
    result, reason = is_thiosugar(test_smiles)
    print(result, "->", reason)