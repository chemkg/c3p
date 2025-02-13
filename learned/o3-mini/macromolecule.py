"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecule 
Definition: “A macromolecule is a molecule of high relative molecular mass, the structure of which essentially comprises the multiple repetition of units derived, actually or conceptually, from molecules of low relative molecular mass.”
Improved Heuristic:
  - The SMILES must be valid.
  - The molecular weight must be ≥ 500 Da and the molecule must have at least 15 heavy atoms.
  - Instead of simply looking at the maximum repeated substructure from a Morgan fingerprint (radius=2),
    we compute a "repetition density": for every fingerprint bit that appears at least 3 times, subtract 2 (the minimum)
    and sum these extras, then divide by the heavy atom count.
  - We further require that at least 2 distinct fingerprint bits occur three or more times.
    
These additions were motivated by the observation that many small or borderline molecules artificially
score “repetition” while genuine macromolecules (such as oligosaccharides and peptides) display many repeated environments.
This heuristic is still approximate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines whether a molecule qualifies as a macromolecule based on a heuristic that
    combines size (molecular weight and heavy atom count) with evidence for repetition of substructural
    motifs derived from a Morgan fingerprint.
    
    Heuristic criteria:
      1. Must parse successfully.
      2. Must have molecular weight ≥ 500 Da and at least 15 heavy atoms.
      3. Compute the Morgan fingerprint (radius 2) and extract counts per bit.
      4. Count the number of distinct bits that appear at least 3 times (repeated environments).
      5. Define "repetition density" = (sum_{for each bit with count ≥ 3}(count – 2)) / (number of heavy atoms).
         This measures the number of “extra” repeats relative to size.
      6. Require that at least 2 distinct bits occur ≥ 3 times and that the repetition density is at least 0.15.
         
    Args:
        smiles (str): SMILES string.
    
    Returns:
        bool: True if the molecule is classified as a macromolecule, False otherwise.
        str: Explanation of the reasoning.
    """
    # Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Compute molecular weight and heavy atom count.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da); macromolecules are generally ≥ 500 Da."
    if heavy_atoms < 15:
        return False, f"Not enough heavy atoms ({heavy_atoms}); molecule likely too small for a macromolecule."
    
    # Compute Morgan fingerprint (radius 2) with counts.
    fp = rdMolDescriptors.GetMorganFingerprint(mol, radius=2)
    counts = fp.GetNonzeroElements()  # dict: bit -> count
    
    # Identify all bits that repeat at least 3 times.
    repeated_bits = {bit: cnt for bit, cnt in counts.items() if cnt >= 3}
    if len(repeated_bits) < 2:
        return False, (f"Insufficient distinct repeating substructural units (only {len(repeated_bits)} found); "
                       "macromolecules usually have multiple repeating units.")
    
    # Compute repetition density: extra repeats beyond a minimum of 2 per bit, normalized by heavy atoms.
    extra_repeats = sum((cnt - 2) for cnt in repeated_bits.values())
    repetition_density = extra_repeats / heavy_atoms
    
    if repetition_density < 0.15:
        return False, (f"Low repetition density ({repetition_density:.2f}); macromolecules usually exhibit a "
                       "higher number of repeated units relative to their size.")
    
    # Passed all tests.
    return True, (f"Molecule qualifies as a macromolecule (Molecular weight: {mol_wt:.2f} Da; "
                  f"Heavy atoms: {heavy_atoms}; distinct repeating bits: {len(repeated_bits)}; "
                  f"repetition density: {repetition_density:.2f}).")

# Example usage (for manual testing):
if __name__ == "__main__":
    test_smiles = [
        # A typical macromolecule example (oligosaccharide/complex repeating unit)
        "O([C@H]1[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H](O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO)O[C@H]4[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]4CO)O)[C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O)CO)[C@@H]5O)CO",
        # A small molecule that should not be flagged as a macromolecule.
        "CC1=CC(=O)NC(=O)C1"
    ]
    for smi in test_smiles:
        flag, reason = is_macromolecule(smi)
        print("SMILES:", smi)
        print("Classification:", flag, ";", reason)
        print("-" * 60)