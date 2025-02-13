"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Monosaccharide
A monosaccharide here is defined as a single (non–glycosidically connected) polyhydroxy aldehyde or ketone 
(with potential to exist as a cyclic hemiacetal), with three or more carbons.
This implementation uses two prongs:
  1. A cyclic sugar test: looks for a 5–membered (furanose) or 6–membered (pyranose) ring that contains exactly one
     ring oxygen and a sufficient number of hydroxyl substituents, plus an overall elemental composition “sugar–like”
     (roughly one oxygen per carbon).
  2. An open–chain sugar test: if no ring is found, the molecule must be small (3–8 carbons, MW <300 Da),
     must show at least one carbonyl group and a high “–OH density” (with slightly relaxed requirements if a carboxyl is present).
If one of these tests passes the molecule is classified as a monosaccharide.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Applies both a cyclic sugar test and (if that fails) an open–chain sugar test.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if classified as a monosaccharide, else False.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, check if the molecule is a single, connected unit.
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        return False, "Molecule contains multiple fragments – likely a glycosidic conjugate."

    # Global composition: count carbons and oxygens.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(carbons) == 0:
        return False, "No carbon atoms found."
    co_ratio = len(oxygens) / len(carbons)
    # For a typical free monosaccharide the ratio is ~1.0 (eg glucose: C6O6). 
    # We require at least 0.9; if too low then extra non‐sugar (or aromatic) carbons are present.
    if co_ratio < 0.9:
        return False, f"Overall O/C ratio is too low ({co_ratio:.2f}); not sugar–like."

    # ----- 1. Cyclic sugar (hemiacetal) test -----
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    for ring in ring_info.AtomRings():
        # Only examine 5–membered (furanose) or 6–membered (pyranose) rings.
        if len(ring) not in (5, 6):
            continue
        oxygen_in_ring = 0
        ring_carbons = []
        for idx in ring:
            a = mol.GetAtomWithIdx(idx)
            if a.GetAtomicNum() == 8:
                oxygen_in_ring += 1
            elif a.GetAtomicNum() == 6:
                ring_carbons.append(idx)
        if oxygen_in_ring != 1:
            continue  # Reject if not exactly one ring oxygen.
        
        # Count --OH substituents on ring carbons.
        OH_count = 0
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            # Check neighbors not in the ring.
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue
                # Using explicit hydrogens (temporarily we add them)
                # A hydroxyl oxygen is an oxygen with at least one hydrogen.
                if nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() >= 1:
                    OH_count += 1
                    break  # Count at most one –OH per ring carbon.
        # For a furanose we require at least 2; for pyranose at least 3.
        if (len(ring) == 5 and OH_count >= 2) or (len(ring) == 6 and OH_count >= 3):
            candidate_rings.append((ring, OH_count))
    
    if len(candidate_rings) > 0:
        # If more than one candidate ring is found, we assume a glycosidic (multi sugar) structure.
        if len(candidate_rings) > 1:
            return False, "Multiple candidate sugar rings detected – likely an oligosaccharide or glycoside."
        ring, oh = candidate_rings[0]
        ring_size = len(ring)
        ring_type = "furanose-like" if ring_size == 5 else "pyranose-like"

        # Additional check: a free monosaccharide (or a single sugar unit) should have a composition close to typical sugars.
        # Many simple sugars (and uronic acids or amino sugars) have relatively few extra heavy atoms. 
        # We check that the total molecular weight is not excessively high.
        mw = rdMolDescriptors.CalcExactMolWt(mol)
        # Typical monosaccharides fall roughly in the 80–300 Da range.
        if not (80 <= mw <= 300):
            return False, f"Identified a {ring_type} ring with {oh} hydroxyl(s) but the overall MW ({mw:.1f} Da) is out of range for a free monosaccharide."
        # Also the overall O/C ratio has already been checked.
        reason = (f"Contains a {ring_type} ring (size {ring_size}) with {oh} hydroxyl substituents on ring carbons "
                  f"and overall MW {mw:.1f} Da (O/C ratio ~{co_ratio:.2f}).")
        return True, reason

    # ----- 2. Open–chain sugar test -----
    # Count carbons: free (open–chain) sugars normally have between 3 and 8 carbons.
    carbon_count = len(carbons)
    if carbon_count < 3:
        return False, f"Too few carbons ({carbon_count}); need at least 3 for a monosaccharide."
    if carbon_count > 8:
        return False, f"Too many carbons ({carbon_count}); open–chain monosaccharides normally have 3–8 carbons."
    
    # Add explicit hydrogens to better count hydroxyls.
    mol_H = Chem.AddHs(mol)
    # Define a SMARTS for hydroxyl group (-OH).
    OH_smarts = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol_H.GetSubstructMatches(OH_smarts)
    hydroxyl_count = len(hydroxyl_matches)
    
    # Check for a carbonyl group (either aldehyde [H]C(=O) or ketone [#6]C(=O)).
    carbonyl_smarts = Chem.MolFromSmarts("[CX3](=O)")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_smarts)
    if not has_carbonyl:
        return False, "No explicit carbonyl detected – cannot classify as an open–chain monosaccharide."
    
    # Use a slightly relaxed molecular weight cutoff for open–chain sugars.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw > 300:
        return False, f"Molecular weight too high for an open–chain monosaccharide ({mw:.1f} Da)."
    
    # Heuristically, many free sugars have almost one hydroxyl per carbon.
    # In an aldose we expect roughly (carbon_count - 1) hydroxyls.
    # (If a carboxyl group is present, one –OH is replaced so we require one less.)
    # We also allow a sugar that is deoxy (losing one –OH) so relax further if necessary.
    # For our purposes, require at least: carbon_count - 2.
    required_oh = carbon_count - 2
    if hydroxyl_count < required_oh:
        return False, f"Only {hydroxyl_count} hydroxyl groups detected but require at least {required_oh} for an open–chain monosaccharide with {carbon_count} carbons."
    
    reason = (f"Molecule is an open–chain monosaccharide with {carbon_count} carbons, "
              f"{hydroxyl_count} hydroxyl groups, MW {mw:.1f} Da, and a carbonyl group present (O/C ratio ~{co_ratio:.2f}).")
    return True, reason

# Uncomment the code below to run some tests:
if __name__ == "__main__":
    test_examples = [
        # True positives:
        ("O1[C@@]([C@H](O)[C@@H](O)[C@@H]1O)([C@H](O)CO)[H]", "beta-D-glucofuranose"),
        ("[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)CO", "aldehydo-L-arabinose"),
        ("[H][C@]1(O[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)CO", "beta-D-galactofuranose"),
        ("OC[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O", "alpha-L-galactose"),
        ("OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@@H](O)[C@@H]1O", "alpha-D-manno-heptulopyranose"),
        ("O1[C@H]([C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O)C", "6-deoxy-beta-L-talopyranose"),
        ("OC[C@@H]1OC(O)[C@@H](O)[C@H]1O", "L-ribofuranose"),
        # Open–chain examples:
        ("OCC(=O)[C@@H](O)[C@H](O)[C@H](O)C([O-])=O", "5-dehydro-D-gluconic acid"),
    ]
    for smi, name in test_examples:
        result, msg = is_monosaccharide(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {msg}\n{'-'*60}")