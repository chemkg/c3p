"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: Alpha-Amino Acid
Definition: An amino acid in which the amino group is located on the carbon atom 
at the position alpha to the carboxy group.

This improved function uses a two‐step strategy:
  1. It identifies a carboxyl group (using a SMARTS pattern that matches C(=O)[O,OH],
     which will capture both protonated and deprotonated free carboxyl groups).
     
  2. For each carboxyl group found, it considers each neighboring carbon as a candidate 
     for the alpha carbon. A valid alpha carbon should have a total valence (explicit 
     heavy atom bonds plus implicit hydrogens) of 4 – this covers standard alpha amino 
     acids (CH with 1 H) as well as glycine (CH2 with 2 H’s). Moreover, one of its neighbors 
     (aside from the carboxyl carbon) must be a nitrogen that does not appear “tied up” in 
     an amide linkage (i.e. not bonded to a second carbonyl group).
     
If these criteria are satisfied then the molecule is classified as an alpha–amino acid.
"""

from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an alpha–amino acid.
        str: A reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a free carboxyl group.
    # This pattern covers both protonated (-OH) and deprotonated (-O) forms.
    pattern_carboxyl = Chem.MolFromSmarts("C(=O)[O,OH]")
    carboxyl_matches = mol.GetSubstructMatches(pattern_carboxyl)
    
    if not carboxyl_matches:
        return False, "No free carboxyl group (C(=O)[O,OH]) found in the molecule"
    
    # For every carboxyl group found, examine its neighboring carbon (candidate alpha carbon)
    for match in carboxyl_matches:
        # In the SMARTS pattern "C(=O)[O,OH]", the first atom (index 0) is the carboxyl carbon.
        carboxyl_idx = match[0]
        carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
        
        # Iterate over neighbors of the carboxyl carbon.
        for neighbor in carboxyl_atom.GetNeighbors():
            # Candidate alpha must be a carbon.
            if neighbor.GetAtomicNum() != 6:
                continue
            
            # Compute the total valence: explicit heavy neighbors + implicit hydrogens.
            total_valence = neighbor.GetDegree() + neighbor.GetTotalNumHs()
            # For most amino acids this equals 4 (CH with 1 H) 
            # while glycine (CH2) gives 2 (neighbors) + 2 implicit H = 4.
            if total_valence != 4:
                continue
            
            # We now want to ensure that one of the other neighbors (besides the carboxyl group)
            # is an amino (nitrogen) substituent that is "free" (i.e. not acylated).
            found_free_amino = False
            for subnbr in neighbor.GetNeighbors():
                # Skip the bond to the carboxyl carbon.
                if subnbr.GetIdx() == carboxyl_idx:
                    continue
                if subnbr.GetAtomicNum() == 7:  # a nitrogen candidate
                    # Check that this nitrogen is not part of an amide (i.e. not bonded to a carbonyl 
                    # carbon apart from the alpha carbon connection). If the nitrogen (aside from the bond 
                    # to the alpha candidate) is attached to any carbon that is double-bonded to oxygen, we
                    # suspect that it is involved in a peptide/amide bond.
                    involved_in_amide = False
                    for nn in subnbr.GetNeighbors():
                        # Skip the connection back to the candidate alpha carbon.
                        if nn.GetIdx() == neighbor.GetIdx():
                            continue
                        if nn.GetAtomicNum() == 6:
                            bond = mol.GetBondBetweenAtoms(subnbr.GetIdx(), nn.GetIdx())
                            if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                                # Found a C=O neighboring the nitrogen (and not the carboxyl group)
                                involved_in_amide = True
                                break
                    if involved_in_amide:
                        continue
                    # If we get here the nitrogen appears free.
                    found_free_amino = True
                    break
            if found_free_amino:
                return True, ("Alpha amino acid pattern found: an alpha carbon (next to a free carboxyl group) "
                              "has a free (non-amide) amino substituent")
    return False, "No alpha amino acid pattern found: no suitable alpha carbon with a free amino group next to a carboxyl group"

# For testing purposes you can run some example SMILES strings.
if __name__ == "__main__":
    test_smiles = [
        "OC(=O)C(N)CN1N=CC=C1",          # 3-(1-Pyrazolyl)-alanine ★ True
        "N[C@@H](CC1=CC=C(F)C=C1)C(O)=O",  # 4-fluorophenyl-L-alanine ★ True
        "CN[C@@H](Cc1ccccc1)C(O)=O",       # N-methyl-L-phenylalanine ★ True
        # An example of a known false‐positive from previous attempts:
        "CNC(=N)NCCC[C@H](N)C([O-])=O",    # N(omega)-methyl-L-argininate (should be False)
        # A few false negatives (for which our criteria should not match):
        "O=C1[C@H]2[C@@](CC[C@H]3[C@H]2[C@@H](O)[C@H](O)[C@H]4[C@@H](O)[C@H](O)CC[C@]34C)(C)[C@H](C1)[C@H]([C@@H]5OC(O)C[C@@H](C5)C(C)C)C",  # Contignasterol: False
    ]
    
    for s in test_smiles:
        result, reason = is_alpha_amino_acid(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")