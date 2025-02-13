"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Monosaccharide
A monosaccharide here is defined as a single (non–glycosidically connected) polyhydroxy aldehyde or ketone 
(with potential to exist as a cyclic hemiacetal), with three or more carbons.
This implementation uses a two–pronged heuristic:
  1. A cyclic sugar test: looks for a furanose (5–membered) or pyranose (6–membered) ring that contains exactly one ring oxygen 
     and a sufficiently high number of hydroxyl (–OH) substituents on the ring carbons.
  2. An open–chain sugar test: if no candidate ring is found then the molecule must be small (3–8 carbons and MW < 250),
     have an explicit carbonyl group and a high density of –OH groups. In addition, if an aldehyde is present we require one more –OH than
     if a ketone is presumed, and we check that most of the oxygens are indeed in hydroxyls.
If one of these tests passes the molecule is classified as a monosaccharide.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Applies a cyclic sugar test (looking for a furanose/pyranose ring having one ring oxygen and enough –OH substitutions)
    and (if that fails) an open–chain sugar test (small molecule with an explicit carbonyl and many hydroxyls).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if classified as a monosaccharide, else False.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure molecule is a single connected unit.
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        return False, "Multiple fragments detected – likely a glycosidic conjugate"
    
    # --- 1. Cyclic sugar test ---
    # Look for candidate rings that are 5– or 6–membered (furanose or pyranose)
    candidate_rings = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6):
            continue  # not the expected size
        
        # Count atoms in the ring: we require exactly one oxygen and the rest carbons.
        oxygen_in_ring = 0
        ring_carbon_idxs = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_in_ring += 1
            elif atom.GetAtomicNum() == 6:
                ring_carbon_idxs.append(idx)
        if oxygen_in_ring != 1:
            continue
        
        # Now count –OH substitutions on ring carbons.
        OH_on_ring = 0
        for idx in ring_carbon_idxs:
            atom = mol.GetAtomWithIdx(idx)
            # Look for an oxygen neighbor (outside the ring) that has at least one hydrogen.
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 8 and nb.GetIdx() not in ring:
                    # Using the molecule with explicit hydrogens to be sure:
                    if nb.GetTotalNumHs() >= 1:
                        OH_on_ring += 1
                        break  # count one –OH per ring carbon
        # For furanoses (5–membered ring) we expect at least 2 ring carbons hydroxylated;
        # for pyranoses (6–membered) at least 3.
        if (len(ring) == 5 and OH_on_ring >= 2) or (len(ring) == 6 and OH_on_ring >= 3):
            candidate_rings.append(ring)
    
    if len(candidate_rings) == 1:
        ring_size = len(candidate_rings[0])
        OH_on_ring = 0
        for idx in candidate_rings[0]:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                for nb in atom.GetNeighbors():
                    if nb.GetAtomicNum() == 8 and nb.GetIdx() not in candidate_rings[0] and nb.GetTotalNumHs() >= 1:
                        OH_on_ring += 1
                        break
        ring_type = "furanose-like" if ring_size == 5 else "pyranose-like"
        reason = (f"Contains a {ring_type} ring (size {ring_size}) with {OH_on_ring} hydroxyl substituents on ring carbons.")
        # If extra sugar rings are found, we assume a glycosidic (multi–sugar) structure.
        if len(candidate_rings) > 1:
            return False, "Multiple candidate sugar rings detected – likely an oligosaccharide/glycoside"
        return True, reason
    
    # --- 2. Open-chain sugar test (when no appropriate cyclic sugar is found) ---
    # Count carbons; a monosaccharide should have at least 3.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, f"Too few carbons ({carbon_count} found, need at least 3)"
    
    # For better –OH detection, add explicit hydrogens.
    mol_with_H = Chem.AddHs(mol)
    OH_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol_with_H.GetSubstructMatches(OH_pattern)
    hydroxyl_count = len(hydroxyl_matches)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Look for an explicit carbonyl group.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    
    # Use a molecular weight cutoff to help exclude larger conjugates.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Open–chain monosaccharides are usually small.
    if mw > 250:
        return False, f"Molecular weight too high for an open–chain monosaccharide ({mw:.1f} Da)"
    
    # Likewise, free sugars typically have no more than 8 carbons.
    if carbon_count > 8:
        return False, f"Too many carbons for an open–chain monosaccharide ({carbon_count} found)"
    
    if not has_carbonyl:
        return False, "No explicit carbonyl detected – cannot assign as open–chain monosaccharide"
    
    # Heuristic on –OH density: many free sugars have nearly one –OH per carbon.
    # Also, if an aldehyde is present (smarts: [H]C(=O)) require slightly more –OH than if a ketone is implied.
    aldehyde_pattern = Chem.MolFromSmarts("[H]C(=O)")
    is_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    req_OH = carbon_count - 1 if is_aldehyde else carbon_count - 2
    if hydroxyl_count < req_OH:
        return False, (f"Not enough hydroxyl groups for an open–chain monosaccharide: "
                       f"{hydroxyl_count} found but require at least {req_OH} (for {carbon_count} carbons)")
    
    # Also check that most oxygens are in –OH groups. (Typical ratios for sugars are high.)
    if oxygen_count > 0 and (hydroxyl_count / oxygen_count) < 0.7:
        return False, ("Low proportion of hydroxyl groups relative to total oxygen count; "
                       "this open–chain molecule is not sugar–like")
    
    reason = (f"Molecule is an open–chain monosaccharide with {carbon_count} carbons, "
              f"{hydroxyl_count} hydroxyl groups, and molecular weight {mw:.1f} Da, "
              f"with an explicit carbonyl group detected.")
    return True, reason

# Example usage (you can uncomment the following lines to test):
# examples = [
#    "O1[C@@]([C@H](O)[C@@H](O)[C@@H]1O)([C@H](O)CO)[H]",  # beta-D-glucofuranose (cyclic)
#    "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)CO",              # aldehydo-L-arabinose (open–chain)
#    "OCC(=O)[C@@H](O)[C@H](O)[C@H](O)C([O-])=O",        # keto-D-fructuronate (false positive in previous attempt)
#    "CC1=C(C=C)\\C(NC1=O)=C/c1[nH]c(Cc2[nH]c(\\C=C3NC(=O)C(C=C)=C\\3C)c(C)c2CCC(O)=O)c(CCC(=O)O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O)c1C",  # mono(glucosyluronic acid)bilirubin (cyclic)
# ]
# for smi in examples:
#     result, msg = is_monosaccharide(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {msg}\n")