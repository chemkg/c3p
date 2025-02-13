"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary Amine
Definition: A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.
For the purposes of this classifier, the molecule is considered a primary amine if it contains at least one –NH2 group 
(where the nitrogen has exactly two hydrogens and one heavy-atom neighbor) that is not part of an amide bond.
Also, if the molecule contains any amide bonds (indicative of peptides or related compounds) it is rejected.
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine (R–NH2) is defined as a derivative of ammonia in which exactly one hydrogen is replaced by a
    hydrocarbyl (or related) substituent. Here we require:
      - The nitrogen has exactly 2 hydrogens.
      - It has exactly 1 non-hydrogen neighbor (the substituent).
      - That neighbor is not part of a carbonyl group (which would indicate an amide bond).
    Also, if the molecule contains amide bonds (typical for peptides), the molecule is not classified as a (simple) primary amine.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a primary amine, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we can count them correctly
    mol = Chem.AddHs(mol)
    
    # First, if the molecule contains any amide bonds it is likely a peptide or amide derivative.
    # We look for a bonded pattern "C(=O)N" (note: this is a heuristic).
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if amide_smarts is None:
        return False, "Error in generating amide SMARTS pattern"
    if mol.HasSubstructMatch(amide_smarts):
        return False, "Molecule contains amide bond(s) typical of peptides/amides"
    
    # Iterate over all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        
        # Count the hydrogens attached to this nitrogen using the explicit model.
        h_count = atom.GetTotalNumHs()
        # List heavy neighbors (non-hydrogen atoms)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        
        # We require exactly 2 hydrogens and exactly 1 heavy neighbor 
        # (i.e. one substituent replacing one hydrogen from ammonia).
        if h_count != 2 or len(heavy_neighbors) != 1:
            continue
        
        # Now, check if the substituent (the heavy neighbor) is acceptable.
        # For a typical primary amine it should be the hydrocarbyl group, usually containing carbon.
        substituent = heavy_neighbors[0]
        # (Optional) Even if the substituent is not a carbon, check if it is attached further to a carbon.
        # For simplicity we allow any substituent, as long as the N is NH2.
        
        # Check for possible amide connectivity at the substituent.
        # If the substituent is carbon, see if it is part of a carbonyl (i.e. double bond to oxygen)
        if substituent.GetAtomicNum() == 6:
            is_amide = False
            for bond in substituent.GetBonds():
                # Make sure the bond connects substituent and oxygen
                # and that the bond is a double bond.
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(substituent)
                    if other.GetAtomicNum() == 8:
                        is_amide = True
                        break
            if is_amide:
                continue  # Skip this candidate since its substituent is part of a carbonyl
        # If we reached here, we consider that this nitrogen qualifies as a primary amine.
        return True, "Molecule contains a valid primary amine (R–NH2) group."
    
    return False, "No valid primary amine group (with two hydrogens and one substituent not involved in an amide) found."

# Example test calls (you may remove these before use)
if __name__ == "__main__":
    test_smiles = [
        "Cc1cc(cc(c1N)S(O)(=O)=O)C(=C1\\C=CC(=N)C(=C1)S(O)(=O)=O)\\c1ccc(N)c(c1)S(O)(=O)=O",  # acid fuchsin (free acid form) [TP]
        "CN",  # methylamine [TP]
        "O=C(N[C@@H](C(C)C)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC=1C=2C(NC1)=CC=CC2)CCCCN"  # a peptide fragment [FP]
    ]
    for s in test_smiles:
        res, reason = is_primary_amine(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")