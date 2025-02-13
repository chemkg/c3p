"""
Classifies: CHEBI:26660 sesterterpenoid
"""
#!/usr/bin/env python3
"""
Classifies: Sesterterpenoid
Definition: Any terpenoid derived from a sesterterpene (original C25 skeleton, possibly rearranged or missing one or more skeletal atoms).
Note: Terpenoid structures may be highly modified. Therefore, our analysis uses several heuristic criteria:
  - Valid SMILES representation.
  - Carbon count: For acyclic molecules we require 23–28 carbons (C25 skeleton); for cyclic ones we allow 20–70.
  - Molecular weight: relaxed to between 250 and 1300 Da.
  - Presence of isoprene-like fragments, detected via SMARTS "C=C(C)"; at least 2 occurrences are required.
  - Aromaticity: a high fraction of aromatic carbons is atypical for a terpenoid; molecules with >30% of their carbons aromatic are rejected.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is likely a sesterterpenoid based on its SMILES string.
    
    Heuristics used:
      - The SMILES must be valid.
      - For acyclic molecules (0 rings), the carbon count must be between 23 and 28; for cyclic molecules (>=1 ring),
        between 20 and 70 carbons.
      - The exact molecular weight must be between 250 and 1300 Da.
      - At least two isoprene-like fragments (SMARTS pattern "C=C(C)") must be present.
      - The fraction of aromatic carbons (out of all carbons) should not exceed 0.30.
      
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if it is likely a sesterterpenoid, False otherwise.
        str: Explanation of the outcome.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    
    # Count aromatic carbon atoms.
    aromatic_c = [atom for atom in c_atoms if atom.GetIsAromatic()]
    arom_frac = (len(aromatic_c) / c_count) if c_count > 0 else 0.0
    if arom_frac > 0.30:
        return False, f"Aromatic carbon fraction {arom_frac:.2f} is too high for a typical aliphatic terpenoid"
    
    # Determine number of rings.
    try:
        ring_count = Chem.GetSSSR(mol)
    except Exception:
        ring_count = 0
    
    # Check carbon count using different windows depending on ring count.
    if ring_count == 0:
        if c_count < 23 or c_count > 28:
            return False, f"Acyclic molecule: carbon count of {c_count} is not within expected range (23–28) for a C25 core"
    else:
        if c_count < 20 or c_count > 70:
            return False, f"Carbon count of {c_count} is not in allowed range (20–70) for a cyclic sesterterpenoid variant"
    
    # Calculate the exact molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 250 or mw > 1300:
        return False, f"Molecular weight of {mw:.1f} Da is outside the expected range (250–1300 Da) for sesterterpenoids"
    
    # Look for isoprene-like fragments using SMARTS "C=C(C)".
    isoprene_pat = Chem.MolFromSmarts("C=C(C)")
    if isoprene_pat is None:
        return False, "Internal error creating SMARTS pattern"
    isoprene_matches = mol.GetSubstructMatches(isoprene_pat)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene-like fragment(s); expected at least 2 for a sesterterpenoid"
    
    return True, "Molecule meets heuristic criteria (carbon count, molecular weight, ring count, isoprene‐like fragments, and low aromaticity) for a sesterterpenoid"

# Example usage; this section can be removed or commented out in production.
if __name__ == "__main__":
    # Test with one of the provided sesterterpenoid SMILES (Ophiobolin I).
    test_smiles = "O=C1C=C(C)[C@@H]2[C@@H]1C(=CC[C@H]3[C@]4(O[C@@H](C=C(C)C)C[C@@H]4C)CC[C@@]3(C2)C)CO"
    result, reason = is_sesterterpenoid(test_smiles)
    print("Is sesterterpenoid:", result)
    print("Reason:", reason)