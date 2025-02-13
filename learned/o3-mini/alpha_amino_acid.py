"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: Alpha-Amino Acid
Definition: An amino acid in which the amino group is located on the carbon atom 
at the position alpha to the carboxy group.
This improved function searches for a free carboxyl group (C(=O)[O;H1,-]) and then:
  1. Identifies a neighboring alpha carbon.
  2. Checks that the alpha carbon has a nitrogen substituent.
  3. Verifies that this nitrogen does NOT appear to participate in a peptide bond 
     (i.e. it is not also attached to an additional carbonyl group that would indicate 
     it is part of an amide linkage between amino acid residues).
If these criteria are met, the molecule is classified as an alpha–amino acid.
"""

from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    The procedure is as follows:
      1. Parse the SMILES string.
      2. Look for a free carboxyl group using the SMARTS "C(=O)[O;H1,-]".
         (This will match carboxy groups not involved in peptide bonds.)
      3. For each free carboxyl group, get its carbon atom (carboxyl carbon) and find
         its carbon neighbor(s) (candidate alpha carbons).
      4. For each alpha candidate, look for a bound amine (nitrogen) substituent.
         Then check that the nitrogen is “free” – meaning that besides the connection 
         to the alpha carbon, it is not attached to another carbonyl carbon (which would 
         indicate an amide/peptide bond).
      5. If such an alpha carbon is found, we classify the molecule as an alpha-amino acid.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an alpha-amino acid, False otherwise.
        str: A reason message for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a free carboxyl group 
    # This matches a carbon double bonded to oxygen and singly bonded to an oxygen that is either protonated or deprotonated.
    pattern_carboxyl = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    carboxyl_matches = mol.GetSubstructMatches(pattern_carboxyl)
    
    if not carboxyl_matches:
        return False, "No free carboxyl group (C(=O)[O,OH]) found in the molecule"
    
    # For every free carboxyl group found
    for match in carboxyl_matches:
        # In the SMARTS "C(=O)[O;H1,-]", the first atom (index 0) is the carboxyl carbon.
        carboxyl_idx = match[0]
        carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
        
        # Find neighboring carbons; these are candidate alpha carbons.
        for neighbor in carboxyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue  # skip if not carbon
            alpha_candidate = neighbor
            # Now, look for a nitrogen substituent on the alpha candidate (the free amine)
            for subnbr in alpha_candidate.GetNeighbors():
                # We skip the carboxyl carbon already considered.
                if subnbr.GetIdx() == carboxyl_idx:
                    continue
                if subnbr.GetAtomicNum() == 7:  # nitrogen found
                    # Check that this nitrogen is not part of an amide linkage.
                    # For each neighbor of the nitrogen (besides the alpha candidate),
                    # if the neighbor is a carbon that is double-bonded to an oxygen (i.e. C=O),
                    # then we suspect a peptide bond.
                    involved_in_amide = False
                    for n_neighbor in subnbr.GetNeighbors():
                        if n_neighbor.GetIdx() == alpha_candidate.GetIdx():
                            continue
                        if n_neighbor.GetAtomicNum() == 6:
                            # Examine bonds of this carbon to see if there is a double bond to oxygen.
                            for bond in n_neighbor.GetBonds():
                                # Check if the bond is a double bond to oxygen.
                                other_atom = bond.GetOtherAtom(n_neighbor)
                                if other_atom.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                                    # We allow the case where n_neighbor is the free carboxyl carbon.
                                    if n_neighbor.GetIdx() != carboxyl_idx:
                                        involved_in_amide = True
                                        break
                            if involved_in_amide:
                                break
                    if involved_in_amide:
                        # The nitrogen appears to be part of a peptide bond.
                        continue
                    
                    # We also check that the nitrogen has at least one hydrogen (i.e. is not fully substituted by acyl groups).
                    has_hydrogen = any(neigh.GetAtomicNum() == 1 for neigh in subnbr.GetNeighbors())
                    if not has_hydrogen:
                        # In some cases (e.g., N-methyl amino acids) the nitrogen may not show a hydrogen,
                        # so we allow that as long as it is not acylated as determined above.
                        pass
                    
                    # If we found an alpha carbon with a bound, free (non-peptidyl) amino group,
                    # we consider this as an alpha amino acid pattern.
                    return True, ("Alpha amino acid pattern found: "
                                  "an alpha carbon (attached to a free carboxyl group) bears an amino substituent "
                                  "that does not appear to be involved in a peptide bond")
    
    # If no alpha carbon candidate with a free amino group was found, then the molecule is not classified as a free alpha-amino acid.
    return False, "No alpha amino acid pattern found: no suitable alpha carbon with a free amino group next to a free carboxyl group"

# For testing purposes, you can run a set of example SMILES strings.
if __name__ == "__main__":
    test_smiles = [
        "OC(=O)C(N)CN1N=CC=C1",          # 3-(1-Pyrazolyl)-alanine (should be True)
        "N[C@@H](CC1=CC=C(F)C=C1)C(O)=O",  # 4-fluorophenyl-L-alanine (should be True)
        "CN[C@@H](Cc1ccccc1)C(O)=O",       # N-methyl-L-phenylalanine (should be True)
        "O=C(N1[C@H](CCC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)CCC(=O)N"  # Asp-Gln-Pro peptide (should be False)
    ]
    
    for s in test_smiles:
        result, reason = is_alpha_amino_acid(s)
        print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")