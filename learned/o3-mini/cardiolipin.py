"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: Cardiolipin 
Definition: A phosphatidylglycerol composed of two molecules of phosphatidic acid covalently linked to a molecule of glycerol.
Heuristic criteria:
  • Exactly 2 phosphorus atoms,
  • At least 4 fatty acid (acyl) ester groups (using a permissive SMARTS pattern),
  • Presence of a glycerol bridging moiety,
  • High molecular weight (typically >1000 Da),
  • Sufficient number of rotatable bonds,
  • And no (or very few) ring systems (cardiolipins are acyclic),
Note: This method is heuristic and may not cover all examples.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    Cardiolipin is defined as a phosphatidylglycerol composed of two phosphatidic acids
    linked via a glycerol molecule.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as cardiolipin, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (1) Check for exactly 2 phosphorus atoms.
    phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if phosphorus_count != 2:
        return False, f"Expected 2 phosphorus atoms, but found {phosphorus_count}"
    
    # (2) Cardiolipins are acyclic. Reject molecules with rings.
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings > 0:
        return False, f"Unexpected ring structure detected (found {num_rings} ring(s)); cardiolipin is acyclic"

    # (3) Check for at least 4 acyl ester groups.
    # Use a more permissive SMARTS pattern that looks for an ester function bound to any carbon.
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 4:
        return False, f"Expected at least 4 acyl ester groups, but found {len(acyl_matches)}"
    
    # (4) Check molecular weight (cardiolipins are large, typically >1000 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low for cardiolipin: {mol_wt:.2f} Da"
    
    # (5) Check that the molecule has a sufficient number of rotatable bonds.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Not enough rotatable bonds ({n_rotatable}) expected for a flexible cardiolipin structure"
    
    # (6) Verify the presence of a central glycerol bridging moiety.
    # We use a simple approximation for a glycerol fragment: HOCH2–CHOH–CH2OH.
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol bridging moiety not found, which is required for cardiolipin"
    
    # If all tests pass, classify the molecule as cardiolipin.
    return True, "Molecule matches heuristic criteria for cardiolipin (2 P atoms, ≥4 acyl chains, glycerol linker, high molecular weight, acyclic)"

# Example usage (when run as a script):
if __name__ == "__main__":
    # Test using one of the provided cardiolipin SMILES examples
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(OC[C@@H](O)COP(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(O)=O)(O)=O"
    result, reason = is_cardiolipin(test_smiles)
    print(f"Result: {result}\nReason: {reason}")