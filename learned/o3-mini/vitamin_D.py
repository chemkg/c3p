"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: Vitamin D compounds (fat-soluble hydroxy seco-steroids)
Definition:
    Any member of a group of fat-soluble hydroxy seco-steroids that exhibit biological
    activity against vitamin D deficiency. Vitamin D can be obtained from sun exposure, food
    and supplements and is biologically inactive and converted into the biologically active 
    calcitriol via double hydroxylation in the body.
    
Heuristic criteria:
 1. Valid molecule from SMILES.
 2. Contains at least one hydroxyl (-OH) group.
 3. Possesses a “seco-steroid” core: typically 3 fused rings rather than the 4 rings found in
    intact steroids.
 4. Has a molecular weight in the typical range for vitamin D compounds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    
    Criteria:
      - The molecule must be valid.
      - It should contain at least one hydroxyl group (-OH).
      - It should have three rings (consistent with a seco-steroid structure).
      - Its molecular weight should be within the expected range.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True and a reason if the molecule fits vitamin D criteria,
                     False and a reason otherwise.
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.1f} Da is out of expected vitamin D range (250-600 Da)"
    
    # Count the number of rings in the molecule
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count != 3:
        return False, f"Expected 3 intact rings for seco-steroid core, found {ring_count}"
    
    # Check for at least one hydroxyl group using a SMARTS pattern for -OH
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl (-OH) group found"
    
    # If the above criteria are satisfied, we assume the molecule is a vitamin D compound.
    return True, "Matches vitamin D criteria: valid secosteroid (3 rings) with hydroxyl groups and proper molecular weight"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Vitamin D3 example (24-fluoro-25-hydroxyvitamin D3)
    test_smiles = "FC(CC[C@H]([C@@H]1[C@@]2([C@@H](CC1)/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)C)C)C(O)(C)C"
    result, reason = is_vitamin_D(test_smiles)
    print(result, reason)