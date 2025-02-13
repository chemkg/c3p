"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecule 
Definition: “A macromolecule is a molecule of high relative molecular mass, the structure of which essentially comprises the multiple repetition of units derived, actually or conceptually, from molecules of low relative molecular mass.”
Heuristic (improved):
   - The SMILES string must parse successfully.
   - Its molecular weight must be at least 500 Da and it must have at least 15 heavy atoms.
   - Its Morgan fingerprint (radius=2) is computed and we require that at least one fingerprint bit appears 3 or more times, 
     which is taken as an indication of repeating subunits.
Note: This is still a heuristic approach.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string) qualifies as a macromolecule.
    Our improved heuristic requires that:
      - The molecule parses successfully.
      - It has a minimal molecular weight (≥ 500 Da) and number of heavy atoms (≥ 15).
      - Its Morgan fingerprint (radius=2) shows evidence of repetition; i.e. some local environment (bit) occurs at least 3 times.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule satisfies our heuristic for being a macromolecule, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Calculate molecular weight and heavy atom count.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    
    # Check minimum size criteria.
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da); macromolecules are generally ≥ 500 Da."
    if heavy_atoms < 15:
        return False, f"Not enough heavy atoms ({heavy_atoms}); molecule likely too small for a macromolecule."
    
    # Compute Morgan fingerprint (with radius 2) with counts.
    fp = rdMolDescriptors.GetMorganFingerprint(mol, radius=2)
    counts = fp.GetNonzeroElements()
    
    # Determine the maximum repetition count of any bit.
    max_repeat = max(counts.values()) if counts else 0
    if max_repeat < 3:
        return False, (f"No sufficiently repeated substructural unit detected (max count = {max_repeat}); "
                       f"macromolecules usually contain repeated units.")
    
    # If all tests are passed, classify as a macromolecule.
    return True, (f"Molecule qualifies as a macromolecule (Molecular weight: {mol_wt:.2f} Da; "
                  f"Heavy atoms: {heavy_atoms}; max repeating unit count: {max_repeat}).")

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        # Example expected to be a macromolecule (e.g. an oligosaccharide-like structure with repetitive units).
        "O([C@H]1[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H](O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO)O[C@H]4[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]4CO)O)[C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O)CO)[C@@H]5O)CO", 
        # Example expected to be non-macromolecular (small molecule with no obvious repetition).
        "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N1)SCCNC(=O)CCNC(=O)C)"
         "C)C"
    ]
    for smi in test_smiles:
        flag, reason = is_macromolecule(smi)
        print("SMILES:", smi)
        print("Classification:", flag, ";", reason)
        print("-" * 60)