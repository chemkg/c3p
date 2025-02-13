"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: Icosanoid
Definition: Any member of the group of signalling molecules arising from oxidation
of the three C20 essential fatty acids (EFAs): icosapentaenoic acid (EPA), arachidonic
acid (AA) and dihomo-gamma-linolenic acid (DGLA).

Heuristic improvements:
  - Parse the SMILES string using RDKit.
  - Check that total number of carbons is in a “reasonable” range (15–60).
  - Require the presence of an oxygenated carbonyl group (as in carboxylic acids or esters).
  - Look for one or more of three core features:
        (a) a five-membered (cyclopentane) ring (determined via ring info),
        (b) a polyene chain consisting of at least three conjugated C=C bonds (using SMARTS),
        (c) the presence of phosphorus (often encountered in phosphorylated derivatives).
  - Also check that the overall number of rings is not excessively high (e.g. > 3) unless a cyclopentane is present.
  - If no phosphate is present, reject molecules with molecular weight > 1000 Da.
If these conditions are met, we classify the molecule as a potential icosanoid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines whether a molecule is an icosanoid based on its SMILES string.
    
    Heuristics:
      - Carbon count must be within 15–60.
      - An oxygenated carbonyl group (e.g. in carboxylic acid or ester) must be present.
      - Moreover, at least one of the following core motifs should be present:
            • a cyclopentane ring (found by examining all rings for one with exactly 5 atoms),
            • a polyene chain with at least three conjugated double bonds,
            • or a phosphorus atom (which often occurs in derivatives such as CDP-DG(…)).
      - Additionally, if no phosphate is present then a very high molecular weight (>1000 Da)
        or an excessive number of rings (> 3) will be used to filter out false positives.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a potential icosanoid, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons (atomic number 6)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (15 <= c_count <= 60):
        return False, f"Carbon count ({c_count}) is not in the expected range (15–60)"
        
    # Check for oxygenated carbonyl group.
    # The SMARTS pattern "[CX3](=O)[OX2,H1]" matches carboxylic acids, esters, etc.
    oxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2,H1]")
    if not mol.HasSubstructMatch(oxy_pattern):
        return False, "No oxygenated carbonyl (e.g. carboxylic acid or ester) found"
    
    # Determine if there is any phosphorus atom (atomic number 15)
    has_phosphate = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    
    # Check for cyclopentane ring by using ring info.
    ring_info = mol.GetRingInfo().AtomRings()
    has_cyclopentane = any(len(ring) == 5 for ring in ring_info)
    
    # Check for a polyene chain: use a SMARTS pattern for at least three conjugated C=C bonds.
    # This pattern "C=C/C=C/C=C" requires three conjugated double bonds.
    polyene_pattern = Chem.MolFromSmarts("C=C/C=C/C=C")
    has_polyene = mol.HasSubstructMatch(polyene_pattern)
    
    # Compute the total number of rings.
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    # For molecules with no phosphate, if the ring count is very high the structure is likely unrelated.
    if not has_phosphate and n_rings > 3:
        return False, f"Too many rings present ({n_rings}) for a typical icosanoid structure"
    
    # Also check molecular weight: if no phosphate is present, very heavy molecules (>1000 Da) are unlikely.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not has_phosphate and mol_wt > 1000:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too high for a typical icosanoid"
    
    # Determine which core motif qualifies.
    core_motif = None
    if has_cyclopentane:
        core_motif = "contains a cyclopentane ring characteristic of prostaglandin cores"
    elif has_polyene:
        core_motif = ("contains a polyene chain indicating several conjugated double bonds "
                      "typical of oxidized fatty acids")
    elif has_phosphate:
        core_motif = ("contains a phosphate group suggesting it is a phosphorylated icosanoid derivative")
    
    if core_motif is None:
        return False, "No cyclopentane ring, sufficient polyene chain, or phosphate group was detected"
    
    # If we reach here, the structure qualifies.
    reason = (f"Carbon count {c_count} within range, oxygenated functionality present, "
              f"molecular weight {mol_wt:.1f} Da, and molecule {core_motif}.")
    return True, reason

# Example usage (for testing):
if __name__ == "__main__":
    # Test with one of the samples from the true positives (Cloprostenol)
    smiles_example = "ClC=1C=C(OC[C@H](O)\\C=C\\[C@@H]2[C@H]([C@@H](O)C[C@H]2O)C/C=C\\CCCC(O)=O)C=CC1"
    result, reason = is_icosanoid(smiles_example)
    print("Is icosanoid:", result)
    print("Reason:", reason)