"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: Icosanoid
Definition: Any member of the group of signalling molecules arising from oxidation
of the three C20 essential fatty acids (EFAs): icosapentaenoic acid (EPA), arachidonic
acid (AA) and dihomo-gamma-linolenic acid (DGLA).

Heuristic improvements:
  - Carbon count between 15 and 60.
  - Require an oxygenated carbonyl group (as in acids/esters).
  - Look for at least one core motif:
       • a cyclopentane ring (using ring information),
       • or a polyene chain of at least three consecutively conjugated double bonds,
       • or a phosphorus atom (commonly seen in phosphorylated derivatives).
  - Relax total ring count restrictions if a cyclopentane ring or phosphate is present;
       otherwise, reject if the number of rings exceeds 5.
  - If no phosphate is present, very heavy molecules (>1000 Da) are unlikely to be icosanoids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines whether a molecule is a potential icosanoid based on its SMILES string.
    
    Heuristics used:
      1. Carbon count must be between 15 and 60.
      2. Molecule must contain an oxygenated carbonyl group (e.g. carboxylic acid or ester).
      3. At least one of these core motifs should be present:
           - A cyclopentane ring (found from ring info),
           - A polyene chain with at least three conjugated C=C bonds (using a modified SMARTS pattern),
           - Or a phosphorus atom (indicative of phosphorylated derivatives).
      4. If neither a cyclopentane ring nor phosphate is present then:
            - the overall ring count should not be excessive (here, > 5 is taken as too many),
            - and the molecular weight must be below 1000 Da.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an icosanoid candidate, False otherwise.
        str: A string detailing the reasoning behind the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons (atomic number 6) and check range 15-60.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (15 <= c_count <= 60):
        return False, f"Carbon count ({c_count}) is not in the expected range (15–60)"
    
    # Check for an oxygenated carbonyl group.
    # SMARTS pattern [CX3](=O)[OX2,H1] should cover carboxylic acids, esters, etc.
    oxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2,H1]")
    if not mol.HasSubstructMatch(oxy_pattern):
        return False, "No oxygenated carbonyl (e.g. carboxylic acid or ester) found"
    
    # Check for a phosphorus atom (atomic number 15)
    has_phosphate = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    
    # Check for cyclopentane ring: Look at ring info for any ring of exactly 5 atoms.
    ring_info = mol.GetRingInfo().AtomRings()
    has_cyclopentane = any(len(ring) == 5 for ring in ring_info)
    
    # Check for a polyene chain: Use a SMARTS pattern for 3 consecutive double bonds.
    # This pattern is written generically (without explicit stereochemical markers)
    # so that molecules with three conjugated C=C bonds are matched.
    polyene_pattern = Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]-[#6]=[#6]")
    has_polyene = mol.HasSubstructMatch(polyene_pattern)
    
    # Compute total number of rings
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    
    # Apply restrictions if neither cyclopentane nor phosphate is present.
    if not (has_cyclopentane or has_phosphate):
        if n_rings > 5:
            return False, f"Too many rings present ({n_rings}) for a typical icosanoid structure without a cyclopentane or phosphate core"
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt > 1000:
            return False, f"Molecular weight ({mol_wt:.1f} Da) is too high for a typical non-phosphorylated icosanoid"
    else:
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Decide on a core motif.
    if has_cyclopentane:
        core_motif = "contains a cyclopentane ring characteristic of prostaglandin cores"
    elif has_polyene:
        core_motif = ("contains a polyene chain indicating several conjugated double bonds typical of oxidized fatty acids")
    elif has_phosphate:
        core_motif = ("contains a phosphate group suggesting it is a phosphorylated icosanoid derivative")
    else:
        # This branch should not be reached because our filtering above would have caught this.
        return False, "No cyclopentane ring, sufficient polyene chain, or phosphate group was detected"
    
    reason = (f"Carbon count {c_count} within range, oxygenated functionality present, " +
              f"molecular weight {mol_wt:.1f} Da, and molecule {core_motif}.")
    return True, reason

# Example usage (for testing):
if __name__ == "__main__":
    # Test using one of the known icosanoid samples: Cloprostenol
    test_smiles = "ClC=1C=C(OC[C@H](O)\\C=C\\[C@@H]2[C@H]([C@@H](O)C[C@H]2O)C/C=C\\CCCC(O)=O)C=CC1"
    result, explanation = is_icosanoid(test_smiles)
    print("Is icosanoid:", result)
    print("Reason:", explanation)