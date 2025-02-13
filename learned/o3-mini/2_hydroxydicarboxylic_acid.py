"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid

Definition: A dicarboxylic acid carrying a hydroxy group on the carbon atom at the position alpha
to (at least one) carboxy group.
Our improved approach: 
 - The molecule must have exactly 2 carboxyl (-C(=O)[OH]) groups.
 - For each carboxyl group, we look at its bonded (non-oxygen) neighbors (“alpha-carbon candidates”).
 - If at least one candidate bears a hydroxy (–OH) substituent then it qualifies.
 - HOWEVER, if all such candidates are sp³ then we require that they be the same carbon (i.e. exactly one unique alpha–OH),
   because many false positives (e.g. tartaric acid) have two separate α–OH centers.
 - If any candidate is sp² we accept that case (to cover cases like dihydroxyfumaric acid).
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
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups show up explicitly
    mol = Chem.AddHs(mol)
    
    # SMARTS for carboxy acid group: C(=O)[OH]
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    # Get unique acid carbon indices (the first atom in the SMARTS)
    acid_carbon_idxs = set(match[0] for match in acid_matches)
    if len(acid_carbon_idxs) != 2:
        return False, f"Expected exactly 2 carboxyl groups; found {len(acid_carbon_idxs)}"
    
    # This set will hold unique α–carbon candidates that carry an OH
    alpha_oh_atoms = set()
    
    # Loop over each acid carbon index.
    for acid_idx in acid_carbon_idxs:
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        # Look at neighbors of the acid carbon; we want carbon neighbors only.
        for neighbor in acid_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:  # ignore non-carbon neighbors
                continue
            # Now check if this neighbor (the potential alpha-carbon) bears an -OH group.
            has_oh = False
            for sub_neigh in neighbor.GetNeighbors():
                if sub_neigh.GetIdx() == acid_idx: 
                    continue  # skip the acid carbon
                if sub_neigh.GetAtomicNum() == 8:
                    # Check if the oxygen is part of an OH (bonded to at least one hydrogen)
                    if any(n.GetAtomicNum() == 1 for n in sub_neigh.GetNeighbors()):
                        has_oh = True
                        break
            if has_oh:
                alpha_oh_atoms.add(neighbor.GetIdx())
    
    if not alpha_oh_atoms:
        return False, "No carboxy group found with an alpha-carbon bearing a hydroxy group"
    
    # For better discrimination:
    # If any of the alpha OH candidates is sp2, we accept the structure.
    for atom_idx in alpha_oh_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetHybridization() == rdchem.HybridizationType.SP2:
            return True, "Found a carboxy group with an alpha-carbon (sp2) bearing a hydroxy group"
    
    # Else, all candidates are sp3.
    # In many true positives (e.g. malic acid derivatives) the same alpha-carbon provides the OH functionality for both COOH groups.
    # If we have more than one unique saturated alpha-carbon, that is similar to the situation in tartaric acid,
    # which we want to exclude.
    if len(alpha_oh_atoms) == 1:
        return True, "Found a carboxy group with an alpha-carbon (sp3) bearing a hydroxy group"
    else:
        return False, f"Found {len(alpha_oh_atoms)} separate saturated alpha-carbons with OH (likely not a 2-hydroxydicarboxylic acid)"

# For testing purposes, one might uncomment and run a few examples:
# test_smiles = [
#     "OC(=O)C(CCC(O)=O)O",           # 2-hydroxyadipic acid (should be True)
#     "O=C(O)CC(O)=O",                # simple dicarboxylic acid but no proper alpha OH (False)
#     "OC(=O)C(\\O)=C(/O)C(O)=O",      # dihydroxyfumaric acid (True, alpha carbon is sp2)
#     "O[C@@H](CC(O)=O)C(O)=O",        # (S)-malic acid (True, one unique sp3 alpha OH)
#     "O[C@@H]([C@H](O)C(O)=O)C(O)=O"   # D-tartaric acid (False: two different alpha-OH centers)
# ]
# for smi in test_smiles:
#     res, reason = is_2_hydroxydicarboxylic_acid(smi)
#     print(f"SMILES: {smi}\nResult: {res}, Reason: {reason}\n")