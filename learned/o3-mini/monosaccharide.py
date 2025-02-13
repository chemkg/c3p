"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: Monosaccharide
A monosaccharide is defined as a single (non‐glycosidically connected) polyhydroxy aldehyde or ketone 
(with potential to exist as a cyclic hemiacetal), having three or more carbon atoms.
This implementation uses a two-pronged heuristic:
  1. A cyclic sugar test: looks for a furanose (5-membered) or pyranose (6-membered) ring that 
     contains exactly one ring oxygen and a sufficient number of ring carbons substituted with hydroxyl groups.
  2. An open-chain sugar test: for non-cyclic forms, the molecule must have an explicit carbonyl group 
     and a high hydroxyl-to-carbon count ratio.
If either condition is met, the function returns True along with the reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    The function applies a cyclic sugar check (furanose/pyranose-like ring with typical –OH substitutions)
    and, if that fails, an open-chain check (explicit carbonyl plus a high density of hydroxyl groups relative 
    to carbon atoms) under the assumption that open-chain sugars have 3 to 8 carbons and many –OH groups.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a monosaccharide, False otherwise
        str: Explanation of the decision
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure molecule is a single fragment (i.e. one connected unit)
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        return False, "Multiple fragments detected – likely a glycosidic conjugate"
    
    # Count carbon atoms; a monosaccharide should have at least 3 carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, f"Too few carbons ({carbon_count} found, need at least 3)"
    
    # To better detect hydroxyl (-OH) groups, add explicit hydrogens.
    mol_with_H = Chem.AddHs(mol)
    # Use SMARTS to locate hydroxyl groups (oxygen with one hydrogen)
    OH_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol_with_H.GetSubstructMatches(OH_pattern)
    hydroxyl_count = len(hydroxyl_matches)
    
    # Look for an explicit carbonyl group (which may indicate an open-chain sugar form)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    
    # -----------------------------------------------------------------
    # 1. Check for cyclic (ring) sugar motif
    # -----------------------------------------------------------------
    ring_found = False
    ring_reason = ""
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6):
            continue  # only consider rings that could be furanose or pyranose
        # Count oxygen atoms in the ring and collect carbon indices
        oxygen_in_ring = 0
        carbon_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_in_ring += 1
            elif atom.GetAtomicNum() == 6:
                carbon_indices.append(idx)
        # For typical sugars, the ring should contain exactly one oxygen.
        if oxygen_in_ring != 1:
            continue
        
        # Now, check for hydroxyl substitutions on ring carbons.
        OH_on_ring = 0
        for idx in carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            # Look through neighbors: if a neighbor is oxygen with an attached hydrogen and is not 
            # part of the ring, consider it an –OH substituent.
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 8 and nb.GetIdx() not in ring:
                    # Count hydrogens on the oxygen; if at least one hydrogen is explicit or implicit.
                    if nb.GetTotalNumHs() >= 1:
                        OH_on_ring += 1
                        break
        # For a 5-membered ring (furanose) there are 4 ring carbons – require at least 2 to be hydroxylated.
        # For a 6-membered ring (pyranose) there are 5 ring carbons – require at least 3 to be hydroxylated.
        if len(ring) == 5 and OH_on_ring >= 2:
            ring_found = True
            ring_reason = f"Contains a furanose-like ring (5-membered) with {OH_on_ring} hydroxyl substitutions on ring carbons."
            break
        elif len(ring) == 6 and OH_on_ring >= 3:
            ring_found = True
            ring_reason = f"Contains a pyranose-like ring (6-membered) with {OH_on_ring} hydroxyl substitutions on ring carbons."
            break
    
    # -----------------------------------------------------------------
    # 2. Check for open-chain sugar motif
    # -----------------------------------------------------------------
    open_chain_found = False
    # We expect open-chain sugars to be relatively low in carbon count (typically 3-8) and heavily hydroxylated.
    if (not ring_found) and has_carbonyl and (carbon_count <= 8):
        # The expected formula for an aldose is roughly C(n)H(2n)O(n), meaning roughly n-1 hydroxyl groups.
        # Allow one less for deoxy sugars.
        if hydroxyl_count >= (carbon_count - 2):
            open_chain_found = True
    
    # -----------------------------------------------------------------
    # Final decision
    # -----------------------------------------------------------------
    if ring_found:
        reason = (f"Molecule is a single unit with {carbon_count} carbons and {hydroxyl_count} hydroxyl groups; " +
                  ring_reason)
        return True, reason
    elif open_chain_found:
        reason = (f"Molecule is an open-chain monosaccharide with {carbon_count} carbons, {hydroxyl_count} hydroxyl groups " +
                  "and an explicit carbonyl group detected.")
        return True, reason
    else:
        return False, "No cyclic sugar ring or open-chain monosaccharide motif detected; not a typical monosaccharide."
        
# Example usage (uncomment to test):
# examples = [
#    "O1[C@@]([C@H](O)[C@@H](O)[C@@H]1O)([C@H](O)CO)[H]",  # beta-D-glucofuranose (cyclic)
#    "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)CO",              # aldehydo-L-arabinose (open-chain)
#    "CCCCC(C)CCCCCCCC\\C=C/CCCCCC[C@@H]1OC[C@H](N)[C@@H]1O", # Penasin E (false negative in previous attempt)
# ]
# for smi in examples:
#    result, msg = is_monosaccharide(smi)
#    print(smi, result, msg)