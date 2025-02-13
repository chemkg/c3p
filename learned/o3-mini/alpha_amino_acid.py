"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: Alpha-Amino Acid
Definition: An amino acid in which the amino group is located on the carbon atom 
at the position alpha to the carboxy group.

This function uses a two‐step strategy:
  1. Identify one or more free carboxyl groups (using a SMARTS pattern that captures both protonated and deprotonated forms).
  2. For each free carboxyl group, examine its neighboring carbon (candidate “alpha carbon”). A valid alpha–carbon is expected to have:
       • A total heavy-atom degree plus implicit H count equal to 4 (covering both CH (e.g. alanine) and CH2 (e.g. glycine)).
       • One additional neighbor (aside from the carboxyl carbon) that is a “free” nitrogen (i.e. not involved in an amide bond).
       
Finally, if exactly one such alpha center is found in the molecule we return True. If none or more than one are found (as in the case of peptides)
we return False with an appropriate reason.
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha–amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a free, monomeric alpha–amino acid.
        str: A reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a free carboxyl group.
    # This should match both the acid (COOH) and deprotonated (COO-) forms.
    pattern_carboxyl = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    carboxyl_matches = mol.GetSubstructMatches(pattern_carboxyl)
    
    if not carboxyl_matches:
        return False, "No free carboxyl group (C(=O)[O,OH]) found in the molecule"
    
    candidate_alpha_count = 0
    candidate_reason = ""
    
    # For each free carboxyl group, look at the carbon it is attached to (candidate alpha carbon)
    for match in carboxyl_matches:
        # In the SMARTS pattern "C(=O)[O;H1,-1]", the first atom (index 0) is the carboxyl carbon.
        carboxyl_idx = match[0]
        carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
        
        # Iterate over neighbors of the carboxyl carbon: possible alpha carbons should be carbons.
        for neigh in carboxyl_atom.GetNeighbors():
            if neigh.GetAtomicNum() != 6:
                continue
            
            # Calculate total valence (explicit heavy neighbors + implicit hydrogens).
            total_valence = neigh.GetDegree() + neigh.GetTotalNumHs()
            if total_valence != 4:
                continue
            
            # Check among the other neighbors of candidate alpha for a "free" amino group.
            found_free_amino = False
            for subnbr in neigh.GetNeighbors():
                # Skip the bond to the carboxyl carbon.
                if subnbr.GetIdx() == carboxyl_idx:
                    continue
                if subnbr.GetAtomicNum() == 7:  # nitrogen candidate
                    # A free amino nitrogen should not be involved in an amide bond.
                    # Here we check that aside from the connection to the candidate alpha carbon,
                    # the nitrogen is not double-bonded to an oxygen.
                    involved_in_amide = False
                    for nn in subnbr.GetNeighbors():
                        if nn.GetIdx() == neigh.GetIdx():
                            continue
                        # If this neighbor of the nitrogen is a carbon that is double-bonded to oxygen, suspect an amide.
                        if nn.GetAtomicNum() == 6:
                            bond = mol.GetBondBetweenAtoms(subnbr.GetIdx(), nn.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                # We suspect an amide involvement.
                                involved_in_amide = True
                                break
                    if involved_in_amide:
                        continue
                    # We found a nitrogen that appears free.
                    found_free_amino = True
                    break  # no need to check further nitrogens
            
            if found_free_amino:
                candidate_alpha_count += 1
                candidate_reason = ("Alpha amino acid pattern found: an alpha carbon (adjacent to a free carboxyl group) "
                                    "has a free (non-amide) amino substituent.")
    
    # If exactly one candidate alpha center is found, we accept the molecule as a free amino acid.
    if candidate_alpha_count == 1:
        return True, candidate_reason
    elif candidate_alpha_count > 1:
        return False, ("Multiple potential alpha-carbon centers were found ({} candidates), "
                       "which suggests the molecule may be a peptide or a derivative rather than a free alpha amino acid."
                       .format(candidate_alpha_count))
    else:
        return False, "No alpha amino acid pattern found: no suitable alpha carbon with a free amino group next to a carboxyl group"

# Simple test runs – you can run this module to see some classifications.
if __name__ == "__main__":
    test_smiles = [
        "OC(=O)C(N)CN1N=CC=C1",          # 3-(1-Pyrazolyl)-alanine (should be True)
        "N[C@@H](CC1=CC=C(F)C=C1)C(O)=O",  # 4-fluorophenyl-L-alanine (should be True)
        "CN[C@@H](Cc1ccccc1)C(O)=O",       # N-methyl-L-phenylalanine (should be True)
        # An example false positive (likely a peptide):
        "O=C(N1[C@H](CCC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)CCC(=O)N",
        # An example false negative (should not match the pattern):
        "C[n+]1ccccc1C([O-])=O"
    ]
    
    for s in test_smiles:
        result, reason = is_alpha_amino_acid(s)
        print(f"SMILES: {s}")
        print(f"Result: {result}")
        print(f"Reason: {reason}\n")