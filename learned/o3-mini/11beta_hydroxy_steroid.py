"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11β-hydroxy steroid
Definition: Any 11-hydroxy steroid in which the hydroxy group at position 11 has beta- configuration.
Note: This is a simplified classifier. Reliable steroid numbering requires a full stereochemical and spatial analysis.
We approximate it by checking that the molecule:
  1. Has a reasonable steroid nucleus (≥17 carbons and several fused rings).
  2. Has at least one beta-configured hydroxy group on a ring carbon (expressed by "[C@@H;R]([O])").
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11β-hydroxy steroid.
    The algorithm first checks that the molecule is steroid-like:
       - at least 17 carbons are present,
       - several rings are present (approximate fused ring system),
    then it attempts to detect a hydroxy group on a ring carbon with beta stereochemistry.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is interpreted as a 11β-hydroxy steroid, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, f"Too few carbon atoms ({carbon_count}) to be a steroid"
    
    # Count rings – steroids usually have 4 fused rings;
    # here we simplify by requiring at least 3 rings.
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 3:
        return False, f"Not enough rings ({ring_count}); typical steroids have 4 fused rings"
    
    # Look for a beta-configured hydroxy group on a ring carbon.
    # The SMARTS "[C@@H;R]([O])" looks for a chiral carbon (with @@ indicating beta configuration in SMILES) 
    # that is a member of a ring (R) and is attached to an -OH group.
    tgt_smarts = "[C@@H;R]([O])"
    tgt_pattern = Chem.MolFromSmarts(tgt_smarts)
    matches = mol.GetSubstructMatches(tgt_pattern)
    
    if not matches:
        return False, "No beta-configured (-O) group found on a ring carbon (candidate for the 11β-hydroxy group)"
    
    # (Optional) Further check the local environment on one candidate match.
    # We expect the carbon to be part of a fused ring: check that it is connected to at least two other ring carbons.
    for match in matches:
        atom = mol.GetAtomWithIdx(match[0])
        # Count how many neighbors are carbons and are in a ring.
        ring_neighbors = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.IsInRing():
                ring_neighbors += 1
        if ring_neighbors >= 2:
            return True, "Molecule has a steroid-like skeleton and contains a beta-configured hydroxy group (candidate for 11β-hydroxy)"
    
    return False, "Found beta-configured hydroxy group candidate(s), but none appear to reside in a fused-ring (steroid) environment"

# Example usage:
if __name__ == "__main__":
    # Test with one of the example SMILES: cortisol.
    test_smiles = "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO"
    result, reason = is_11beta_hydroxy_steroid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")