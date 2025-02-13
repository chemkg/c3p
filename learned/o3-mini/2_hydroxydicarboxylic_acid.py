"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid

Definition: A dicarboxylic acid carrying a hydroxy group on the carbon atom at position alpha
to at least one carboxy group. Our improved strategy is to insist that:
 - The molecule has exactly 2 carboxylic acid groups (–C(=O)OH).
 - For each of these groups, we obtain the carbon neighbor(s) (the potential α–carbon) that bear an –OH group.
 - If any candidate is sp2 we accept (this covers cases such as dihydroxyfumaric acid).
 - Otherwise (if candidates are all sp3) we require that the same carbon is the α–carbon for both acid groups.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 2-hydroxydicarboxylic acid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES and add explicit hydrogens so that -OH groups appear explicitly.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a carboxylic acid group: [CX3](=O)[OX2H]
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    # Get unique indices corresponding to the carbon of the carboxyl groups (first atom in SMARTS)
    acid_carbon_idxs = {match[0] for match in acid_matches}
    if len(acid_carbon_idxs) != 2:
        return False, f"Expected exactly 2 carboxyl groups; found {len(acid_carbon_idxs)}"
    
    # For each carboxyl carbon, find neighboring carbon atoms that carry at least one -OH.
    # For a neighbor to be a valid alpha-carbon candidate it must:
    #  - be a carbon (atomic number 6)
    #  - have at least one neighbor (other than the acid carbon) which is oxygen bound to at least one hydrogen
    def get_alpha_candidates(acid_idx):
        candidates = set()
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        for neighbor in acid_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:  # only consider carbon neighbors
                continue
            # Examine the neighbors of this candidate (excluding the acid carbon)
            for sub_neigh in neighbor.GetNeighbors():
                if sub_neigh.GetIdx() == acid_idx:
                    continue
                if sub_neigh.GetAtomicNum() == 8:  # candidate oxygen
                    # Check if this oxygen has at least one hydrogen attached (i.e. -OH)
                    if any(n.GetAtomicNum() == 1 for n in sub_neigh.GetNeighbors()):
                        candidates.add(neighbor.GetIdx())
                        break
        return candidates

    # Create a mapping from each acid carbon to its set of candidate alpha-carbons
    acid_to_candidates = {}
    for acid_idx in acid_carbon_idxs:
        acid_to_candidates[acid_idx] = get_alpha_candidates(acid_idx)
    
    # Check if at least one acid group has a candidate; if none, then there's no alpha-OH.
    if all(len(cands) == 0 for cands in acid_to_candidates.values()):
        return False, "No carboxy group found with an alpha-carbon bearing a hydroxy group"
    
    # Gather candidates per acid
    candidate_sets = list(acid_to_candidates.values())
    # Compute the union and intersection of candidate sets.
    union_candidates = set().union(*candidate_sets)
    common_candidates = set(candidate_sets[0]).intersection(*candidate_sets[1:])
    
    # If any candidate in the union is sp2 hybridized, accept the structure (covers unsaturated cases).
    for cand_idx in union_candidates:
        cand_atom = mol.GetAtomWithIdx(cand_idx)
        if cand_atom.GetHybridization() == rdchem.HybridizationType.SP2:
            return True, "Found a carboxy group with an alpha-carbon (sp2) bearing a hydroxy group"
    
    # If all candidates are sp3, require that both acid groups share exactly one common candidate.
    if len(common_candidates) == 1:
        return True, "Found a carboxy group with a common alpha-carbon (sp3) bearing a hydroxy group"
    elif len(common_candidates) > 1:
        return False, f"Found {len(common_candidates)} common saturated alpha-carbons with OH (likely not a 2-hydroxydicarboxylic acid)"
    else:
        # There is no common candidate among the two acid groups; hence the alpha-OH centers are separate.
        return False, f"Found separate saturated alpha-carbons with OH for each carboxyl group, not a single common center"

# Example usage (for testing):
# test_smiles = [
#     "OC(=O)C(\\O)=C(/O)C(O)=O",      # dihydroxyfumaric acid: sp2 case, should be True
#     "OC(CCCC(O)=O)C(O)=O",           # 2-hydroxyadipic acid: common sp3 candidate, should be True
#     "O[C@@H](CC(O)=O)C(O)=O",         # (S)-malic acid: common sp3 candidate, should be True
#     "O[C@@H]([C@H](O)C(O)=O)C(O)=O",   # tartaric acid: distinct saturated alpha-OH centers, should be False
# ]
# for smi in test_smiles:
#     result, reason = is_2_hydroxydicarboxylic_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")