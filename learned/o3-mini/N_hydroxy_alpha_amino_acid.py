"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.
Example structures include N-hydroxyglycine, N(6)-hydroxy-L-lysine, N,N-dihydroxy-L-isoleucine, and others.
"""

from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule (provided as a SMILES string) is an N-hydroxy-alpha-amino-acid.
    
    For this classification, we require:
      1. Presence of an amino acid backbone (an alpha-carbon connected to a carboxyl group)
      2. That for at least one amino nitrogen in the backbone, one of its substituents (other than the alpha carbon) is an -OH group
         (i.e. the nitrogen bears a hydroxy substituent).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an N-hydroxy-alpha-amino-acid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Helper function: Is a carbon part of a carboxyl group?
    def has_carboxyl(carbon):
        # A carboxyl carbon typically has a double bond to an oxygen (C=O) and a single bond to another oxygen (O-H).
        doubleO = None
        singleO = None
        for bond in carbon.GetBonds():
            nbr = bond.GetOtherAtom(carbon)
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    doubleO = nbr
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    singleO = nbr
        return (doubleO is not None) and (singleO is not None)
    
    # Now scan through all nitrogen atoms in the molecule.
    # We want to find an amino nitrogen that is part of an amino acid backbone.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # not a nitrogen
        # For each nitrogen, check its neighbors for at least one oxygen (for N-hydroxy)
        n_hydroxy_found = False
        # Also, determine if this nitrogen is in an amino acid backbone by checking for an adjacent alpha carbon.
        backbone_found = False
        alpha_carbon = None
        
        for nbr in atom.GetNeighbors():
            # Look for a carbon candidate (possible alpha carbon)
            if nbr.GetAtomicNum() == 6:
                # Check if this carbon looks like an alpha carbon by having a carboxyl group.
                if has_carboxyl(nbr):
                    backbone_found = True
                    alpha_carbon = nbr
                    break
        # If no backbone connection is detected, skip this nitrogen.
        if not backbone_found:
            continue
        
        # Now check if this nitrogen has an -OH substituent.
        for nbr in atom.GetNeighbors():
            # Exclude the alpha carbon from consideration, because that bond is part of the backbone.
            if alpha_carbon is not None and nbr.GetIdx() == alpha_carbon.GetIdx():
                continue
            # Look for oxygen with a single bond.
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Check if the oxygen likely is part of a hydroxyl (has at least one attached hydrogen).
                # (Implicit H count is used so that even if hydrogens are not explicitly drawn, they are counted.)
                if nbr.GetTotalNumHs() > 0:
                    n_hydroxy_found = True
                    break
                    
        # If we found a backbone nitrogen that carries an -OH substituent, classify as positive.
        if backbone_found and n_hydroxy_found:
            return True, "Amino acid backbone detected and an N-hydroxy substituent found on the amino nitrogen."
    
    return False, "Either no amino acid backbone was detected or no N-hydroxy substitution was found on any amino nitrogen."

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "C(N)(=NO)NCCC[C@H](N)C(=O)O",   # N(5)-[amino(hydroxyimino)methyl]-L-ornithine
        "CSCCCCC(N(O)O)C(O)=O",          # N,N-dihydroxydihomomethionine
        "O=C(O)[C@@H](NO)CCCCCSC",       # N-hydroxy-L-trihomomethionine
        "N[C@@H](CCCCNO)C(O)=O",         # N(6)-hydroxy-L-lysine
        "CC(C)[C@H](N(O)O)C(O)=O",        # N,N-dihydroxy-L-isoleucine
        "ONCC(O)=O",                    # N-hydroxyglycine
        "N[C@@H](CCCC)C(O)=O"            # Amino acid lacking N-hydroxy substitution
    ]
    for smi in test_smiles:
        result, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")