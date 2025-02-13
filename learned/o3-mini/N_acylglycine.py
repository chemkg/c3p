"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N-acylglycine
Definition: An N-acyl-amino acid in which the amino acid specified is glycine.
The key motif is: R-C(=O)-N-CH2-C(=O)O.
This implementation iterates over amide nitrogens to identify a terminal amide
bond where one neighbor is a carbonyl (the acyl group) and the other neighbor is a CH2
(which is further connected to a carboxylic acid group).
"""

from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    The molecule must contain an acyl group attached to a glycine unit via an amide bond,
    i.e., the key motif: R-C(=O)-N-CH2-C(=O)O.
    
    The function identifies candidate amide nitrogens that are terminal (bonded to exactly
    2 heavy atoms) and then tests if one substituent is a carbonyl carbon (acyl group)
    and the other is a CH2 group that is connected to a carboxylic acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a valid N-acylglycine motif is detected, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for proper hydrogen counting.
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms looking for nitrogen atoms that could form an amide.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Skip non-nitrogen atoms.
            
        # For our candidate amide nitrogen, check that it is terminal
        # (i.e. only two heavy-atom neighbors).
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 2:
            continue  # Skip if nitrogen is substituted further (e.g., peptide bonds).
        
        # Consider the two heavy neighbors in both possible assignments:
        # one must be the acyl carbon (part of R-C(=O)) and the other the glycine alpha carbon.
        for acyl_candidate, alpha_candidate in [(heavy_neighbors[0], heavy_neighbors[1]),
                                                (heavy_neighbors[1], heavy_neighbors[0])]:
            # Check acyl candidate: It should be a carbon and look like a carbonyl.
            if acyl_candidate.GetAtomicNum() != 6:
                continue
            # Look among its neighbors (except the nitrogen) for at least one oxygen
            # with a double bond.
            found_carbonyl = False
            for nbr in acyl_candidate.GetNeighbors():
                # Exclude the amide nitrogen.
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                # Check that neighbor is oxygen.
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(acyl_candidate.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                        found_carbonyl = True
                        break
            if not found_carbonyl:
                continue  # Not a valid acyl carbon.
            
            # Check that alpha candidate is a carbon and has exactly 2 hydrogens (CH2).
            if alpha_candidate.GetAtomicNum() != 6:
                continue
            if alpha_candidate.GetTotalNumHs() != 2:
                continue  # Not glycine if not CH2.
            
            # Now check for the glycine carboxyl terminus:
            # The alpha carbon should be bonded to a carboxyl carbon (other than the amide nitrogen).
            acid_candidate = None
            for nbr in alpha_candidate.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue  # Skip amide nitrogen.
                if nbr.GetAtomicNum() == 6:
                    acid_candidate = nbr
                    break
            if acid_candidate is None:
                continue  # No carboxylic acid carbon found.
            
            # Examine the candidate acid carbon.
            # It should be bound to alpha_candidate (already true) and two oxygens: one double bonded and one single bonded to an H.
            oxy_double = 0
            oxy_single = 0
            for nbr in acid_candidate.GetNeighbors():
                # Skip the connection back to alpha carbon.
                if nbr.GetIdx() == alpha_candidate.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(acid_candidate.GetIdx(), nbr.GetIdx())
                    if bond is None:
                        continue
                    # Count based on bond order.
                    if bond.GetBondTypeAsDouble() == 2.0:
                        oxy_double += 1
                    elif bond.GetBondType().name == "SINGLE":
                        # Ensure that this oxygen has at least one hydrogen in the explicit representation.
                        if nbr.GetTotalNumHs() >= 1:
                            oxy_single += 1
            if (oxy_double >= 1 and oxy_single >= 1):
                # We found a valid candidate.
                return True, "Molecule contains the N-acylglycine motif (R-C(=O)-N-CH2-C(=O)O)."
    
    return False, "The molecule does not contain a valid N-acylglycine moiety after filtering."

# Example test cases (you may add more tests if desired).
if __name__ == "__main__":
    test_examples = [
        ("Cepaciachelin", "O=C(N[C@H](C(=O)NCCCCN)CCCCNC(=O)C1=C(O)C(O)=CC=C1)C2=C(O)C(O)=CC=C2"),
        ("N-benzoylglycine", "OC(=O)CNC(=O)c1ccccc1"),
        ("N-(2-furoyl)glycine", "O1C=CC=C1C(NCC(O)=O)=O"),
        ("N-cinnamoylglycine", "OC(=O)CNC(=O)\\C=C\\c1ccccc1"),
        ("butyrylglycine", "C(CNC(=O)CCC)(=O)O"),
        ("N-acetylglycine", "CC(=O)NCC(O)=O"),
    ]
    
    for name, smi in test_examples:
        res, reason = is_N_acylglycine(smi)
        print(f"Name: {name}\nSMILES: {smi}\nClassification: {res}\nReason: {reason}\n")