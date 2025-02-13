"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: CHEBI:2-hydroxydicarboxylic acid

Definition: A dicarboxylic acid carrying a hydroxy group on the carbon atom 
at position alpha to a carboxy group. Our strategy:
 1. The molecule must have exactly 2 carboxyl (–C(=O)OH) groups.
 2. For each carboxyl group, we look at its non-acid neighbors that are carbons.
    A neighbor is accepted as an alpha–carbon candidate if it carries at least one -OH 
    group (detected by an oxygen neighbor that in turn is attached to a hydrogen).
 3. If for any acid group the candidate is unsaturated (sp2), then we require that 
    the same candidate (by index) appears for both carboxyl groups.
 4. Otherwise (if candidates are saturated, sp3), we require that each acid group 
    has at least one candidate – even if they differ.
 
This strategy tries to avoid false positives when unsaturated (sp2) candidates are found 
only on one acid group (as seen in previous attempts), while still accepting many malic acid 
derivatives when each acid “sees” a saturated alpha–OH.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid.
    
    The algorithm:
      - Parses the SMILES string and adds explicit hydrogens.
      - Uses a SMARTS pattern to locate carboxylic acid groups ([CX3](=O)[OX2H]).
      - For each acid group, finds neighboring carbons that have at least one -OH group.
      - If an alpha candidate is sp2, then that candidate must be common to both acids.
        Otherwise, if the candidates are sp3, acceptance is based on each acid having a candidate.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 2-hydroxydicarboxylic acid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES and add explicit hydrogens so that -OH groups are visible.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a carboxylic acid group: [CX3](=O)[OX2H]
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    # Get the unique acid carbon indices (the first atom in the SMARTS match)
    acid_carbon_idxs = {match[0] for match in acid_matches}
    if len(acid_carbon_idxs) != 2:
        return False, f"Expected exactly 2 carboxyl groups; found {len(acid_carbon_idxs)}"
    
    # For each acid carbon, find neighboring carbon atoms that carry at least one -OH group.
    # A neighbor is considered only if it is carbon and (excluding the acid oxygen)
    def get_alpha_candidates(acid_idx):
        candidates = set()
        for neighbor in mol.GetAtomWithIdx(acid_idx).GetNeighbors():
            # Only consider carbon neighbors
            if neighbor.GetAtomicNum() != 6:
                continue
            # Look among the neighbors of this candidate (excluding the acid carbon itself)
            for sub_neigh in neighbor.GetNeighbors():
                if sub_neigh.GetIdx() == acid_idx:
                    continue
                # Check if this neighbor is oxygen and has at least one hydrogen (i.e. -OH)
                if sub_neigh.GetAtomicNum() == 8:
                    if any(n.GetAtomicNum() == 1 for n in sub_neigh.GetNeighbors()):
                        candidates.add(neighbor.GetIdx())
                        break
        return candidates

    # Build a mapping: acid carbon index -> set of candidate alpha-carbon indices.
    acid_to_candidates = {}
    for acid_idx in acid_carbon_idxs:
        acid_to_candidates[acid_idx] = get_alpha_candidates(acid_idx)
    
    # If for any carboxyl group no alpha candidate is found, reject.
    for acid_idx, cands in acid_to_candidates.items():
        if not cands:
            return False, f"No alpha–carbon with -OH found adjacent to acid group at index {acid_idx}"
    
    # For each acid group, partition candidates by hybridization type (sp2 vs sp3).
    acid_to_sp2 = {}
    acid_to_sp3 = {}
    for acid_idx, cands in acid_to_candidates.items():
        sp2_set = set()
        sp3_set = set()
        for cand_idx in cands:
            cand_atom = mol.GetAtomWithIdx(cand_idx)
            if cand_atom.GetHybridization() == rdchem.HybridizationType.SP2:
                sp2_set.add(cand_idx)
            elif cand_atom.GetHybridization() == rdchem.HybridizationType.SP3:
                sp3_set.add(cand_idx)
            else:
                # For other hybridizations we add them to sp3 bucket by default.
                sp3_set.add(cand_idx)
        acid_to_sp2[acid_idx] = sp2_set
        acid_to_sp3[acid_idx] = sp3_set

    # For decision, collect the candidate sets across the two acid groups.
    candidate_sets = list(acid_to_candidates.values())
    common_candidates = set.intersection(*candidate_sets)
    # Also gather the sp2 candidates common to both acids.
    sp2_sets = list(acid_to_sp2.values())
    common_sp2 = set.intersection(*sp2_sets) if sp2_sets else set()
    
    # Now apply decision rules:
    # Rule 1: If any acid group provides only sp2 candidate(s) then require a common sp2 candidate.
    only_sp2 = []
    for acid_idx in acid_carbon_idxs:
        if acid_to_candidates[acid_idx] and (not acid_to_sp3[acid_idx]): 
            only_sp2.append(acid_idx)
    if only_sp2:
        if common_sp2:
            return True, "Found a common unsaturated (sp2) alpha–carbon bearing -OH for the acid group(s)"
        else:
            return False, "Acid groups with unsaturated alpha–carbons do not share a common candidate"
    
    # Rule 2: If there is any common candidate (regardless of hybridization), accept.
    if common_candidates:
        return True, "Found a common alpha–carbon bearing -OH for both carboxyl groups"
    
    # Rule 3: Otherwise, if each carboxyl group has at least one (saturated) alpha–carbon candidate (even if different),
    # we accept the molecule as a 2-hydroxydicarboxylic acid.
    for acid_idx, cands in acid_to_candidates.items():
        if not cands:
            return False, f"No alpha–carbon with -OH found for acid group at index {acid_idx}"
    return True, "Found separate saturated alpha–carbon candidates with -OH for each carboxyl group"

# Example usage (for testing):
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        ("OC(=O)C(\\O)=C(/O)C(O)=O", "dihydroxyfumaric acid"),
        ("OC(=O)\\C=C\\C=C(/O)C(O)=O", "(2Z,4E)-2-hydroxymuconic acid"),
        ("OC(=O)/C=C(/O)C(O)=O",    "(2E,4Z)-2-hydroxymuconic acid"),
        # False positives:
        ("[H]C(=O)C(\\O)=C/C=C(\\O)C(O)=O", "5-formyl-2-hydroxyhepta-2,4-dienedioic acid"),
        ("OC(=O)C(\\O)=C/C(=C\\C=O)C(O)=O", "4-carboxy-2-hydroxy-cis,cis-muconate 6-semialdehyde"),
        # False negatives (should be accepted in our ideal class definition):
        ("OC(CCCC(O)=O)C(O)=O", "2-hydroxyadipic acid"),
        ("O[C@@H](CCC(O)=O)C(O)=O", "(S)-2-hydroxyglutaric acid"),
    ]
    for smi, name in test_smiles:
        result, reason = is_2_hydroxydicarboxylic_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result}, Reason: {reason}\n")