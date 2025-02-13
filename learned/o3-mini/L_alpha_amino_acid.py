"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
#!/usr/bin/env python3
"""
Classifies: L-alpha-amino acid
Definition: Any alpha-amino acid having L-configuration at the alpha-carbon.
An L-alpha-amino acid must have an amino group attached directly to a chiral alpha carbon,
which is also bonded to a carboxylic acid group. For most amino acids (except cysteine),
the natural L form corresponds to an S configuration (by CIP rules). For cysteine the priorities
flip and the natural L form is R.
"""

from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha amino acid based on its SMILES string.
    The algorithm:
      1. Parse the SMILES into an RDKit Mol.
      2. Search for the typical alpha-amino acid motif, 
         i.e., an amino group (N) bonded to a chiral alpha-carbon which in turn is bonded to a carboxyl group.
         Here we attempt to match either "N[C@H](*)C(=O)O" or "N[C@@H](*)C(=O)O". 
      3. For the matched candidate alpha-carbon, use RDKit to assign the stereochemistry.
      4. Identify (heuristically) if the alpha-carbon’s side chain is cysteine-like (i.e. contains sulfur nearby).
         If so, the expected configuration for L is “R”. Otherwise (for most amino acids) we expect “S”.
      5. Return True if the found configuration matches the expected configuration; otherwise, return False with a reason.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where the boolean is True if the molecule matches the class
                     (L-alpha-amino acid) and the reason indicates the outcome.
                     If the SMILES is invalid or motif not found, returns (False, reason).
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns to capture the common motif.
    # The pattern looks for an amino group bonded to a chiral carbon bonded to a carboxylic acid.
    pattern1 = Chem.MolFromSmarts("N[C@H](*)C(=O)O")
    pattern2 = Chem.MolFromSmarts("N[C@@H](*)C(=O)O")
    
    # Try to match either pattern
    matches = mol.GetSubstructMatches(pattern1)
    if not matches:
        matches = mol.GetSubstructMatches(pattern2)
    if not matches:
        return False, "No alpha-amino acid motif found (expected pattern: N[C@H](*)C(=O)O or N[C@@H](*)C(=O)O)"
    
    # Use the first match found.
    # In our SMARTS the order is: index0 = amino N, index1 = chiral alpha carbon, index2 = carboxyl carbon.
    alpha_idx = matches[0][1]
    
    # Now, compute the chiral centers on the molecule.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    # Look for our candidate alpha carbon in the computed chiral centers.
    alpha_center_data = None
    for idx, config in chiral_centers:
        if idx == alpha_idx:
            alpha_center_data = (idx, config)
            break
    if alpha_center_data is None:
        return False, "The candidate alpha-carbon does not have an assigned chiral configuration"
    
    actual_config = alpha_center_data[1]  # 'R' or 'S'
    
    # Identify the side chain: for the alpha-carbon, we expect three heavy-atom neighbors:
    # one is the amino nitrogen (from the pattern) and one is the carboxyl carbon.
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    neighbor_idxs = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors()]
    
    # Remove the known amino nitrogen and carboxyl carbon from neighbor list.
    # Identify the amino neighbor from our match: pattern match index 0.
    amino_idx = matches[0][0]
    # Identify the carboxyl neighbor from our match: pattern match index 2.
    carboxyl_idx = matches[0][2]
    
    side_chain_candidates = [idx for idx in neighbor_idxs if idx not in (amino_idx, carboxyl_idx)]
    
    # If there is no side chain candidate, then the structure is likely glycine (which is achiral) or an error.
    if not side_chain_candidates:
        return False, "No side chain found; glycine is achiral and thus not classified as L-alpha-amino acid"
    
    side_chain_idx = side_chain_candidates[0]
    side_chain_atom = mol.GetAtomWithIdx(side_chain_idx)
    
    # Heuristically determine if the side chain is cysteine-like.
    # In cysteine, the side chain is CH2-SH so the immediate neighbor (a carbon) will be bound to a sulfur atom.
    is_cysteine = False
    if side_chain_atom.GetAtomicNum() == 6:  # carbon
        # Look over the neighbors of the side chain atom (except the alpha-carbon)
        for nbr in side_chain_atom.GetNeighbors():
            if nbr.GetIdx() == alpha_idx:
                continue
            if nbr.GetAtomicNum() == 16:  # sulfur
                is_cysteine = True
                break

    # Determine expected configuration:
    # For most amino acids, L corresponds to an S configuration.
    # For cysteine, due to its sulfur substitution the expected absolute configuration is R.
    expected_config = 'R' if is_cysteine else 'S'
    
    if actual_config == expected_config:
        return True, "Alpha-amino acid with the correct L-configuration found"
    else:
        return (False, f"Alpha-carbon configuration is {actual_config}, but expected {expected_config} for an L-alpha-amino acid")
        
# Example (uncomment to test):
#if __name__ == '__main__':
#    # Test on L-valine SMILES
#    res, reason = is_L_alpha_amino_acid("CC(C)[C@H](N)C(O)=O")
#    print(res, reason)