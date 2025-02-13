"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.
Amino acids normally have an alpha carbon (a carbon attached to exactly three heavy atoms) with one neighbor 
being an amino nitrogen and one neighbor being the carboxyl carbon (which shows a characteristic C(=O)[OH] motif).
The amino nitrogen then must contain at least one hydroxyl substituent.
"""

from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule (provided as a SMILES string) is an N-hydroxy-alpha-amino-acid.
    
    For a molecule to qualify we require:
      1. The presence of an amino acid backbone:
         - An alpha carbon (sp3 carbon) having exactly three heavy-atom neighbors:
           one is an amino nitrogen, one is a carboxyl carbon, and one is the side chain.
         - The carboxyl carbon is defined as a carbon that is doubly bonded to an oxygen
           and singly bonded to an -OH group.
      2. For the amino nitrogen (the one attached to the α-carbon), at least one substituent 
         other than the alpha carbon is an -OH group (i.e. an oxygen with a single bond and at least one hydrogen).
    
    If both conditions are satisfied the function returns True along with an explanation.
    Otherwise it returns False and a reason.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is an N-hydroxy-alpha-amino-acid, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Add explicit hydrogens to help detect -OH (and the H count on oxygens)
    mol = Chem.AddHs(mol)
    
    # Helper function to test if a carbon atom represents a carboxyl carbon.
    def is_carboxyl_carbon(carbon):
        # carboxyl carbon should be sp2 and connected to at least two oxygens,
        # one by a double bond and one by a single bond from an -OH group.
        if carbon.GetAtomicNum() != 6:
            return False
        doubleO = False
        singleOH = False
        for nbr in carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                # detect C=O
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    doubleO = True
                # detect O-H where oxygen is connected by a single bond
                elif bond.GetBondType() == Chem.BondType.SINGLE and nbr.GetTotalNumHs() > 0:
                    singleOH = True
        return doubleO and singleOH

    # Helper function to detect if the amino nitrogen has at least one -OH substituent.
    def has_N_hydroxy(nitrogen, alpha_idx):
        # Check each neighbor except the one from the alpha carbon.
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetIdx() == alpha_idx:
                continue
            if nbr.GetAtomicNum() == 8:  # oxygen candidate
                bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # if oxygen has any hydrogen attached it is likely an -OH group.
                if nbr.GetTotalNumHs() > 0:
                    return True
        return False

    # Iterate over carbons to find a candidate alpha carbon.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # look only at carbons
            continue
        
        # Identify heavy-atom (atomic number > 1) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # In a typical amino acid, the alpha carbon has exactly three heavy neighbors:
        # one amino nitrogen, one carboxyl carbon, and one side-chain atom.
        if len(heavy_neighbors) != 3:
            continue
        
        # Look for an amino nitrogen and a carboxyl carbon among the neighbors.
        amino_nitrogen = None
        carboxyl_carbon = None
        for nbr in heavy_neighbors:
            if nbr.GetAtomicNum() == 7 and amino_nitrogen is None:
                amino_nitrogen = nbr
            elif nbr.GetAtomicNum() == 6:
                # Check if this neighbor qualifies as a carboxyl carbon.
                if is_carboxyl_carbon(nbr):
                    carboxyl_carbon = nbr
        if amino_nitrogen is None or carboxyl_carbon is None:
            continue  # this carbon does not appear to be the alpha carbon of an amino acid
        
        # Check that the nitrogen attached to the α-carbon has an -OH substituent (not counting the bond to the α-carbon).
        if has_N_hydroxy(amino_nitrogen, atom.GetIdx()):
            return True, ("Amino acid backbone detected (α-carbon with exactly three heavy neighbors including a carboxyl group) "
                          "and the amino nitrogen carries an -OH substituent.")
    
    # If we reach here, no candidate alpha carbon with an N-hydroxy amino group was found.
    return False, "No valid N-hydroxy-alpha-amino-acid backbone with the requisite N-hydroxy substitution detected."


# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "C(N)(=NO)NCCC[C@H](N)C(=O)O",    # N(5)-[amino(hydroxyimino)methyl]-L-ornithine (should be positive)
        "CSCCCCC(N(O)O)C(O)=O",           # N,N-dihydroxydihomomethionine (should be positive)
        "O=C(O)[C@@H](NO)CCCCCSC",        # N-hydroxy-L-trihomomethionine (should be positive)
        "N[C@@H](CCCCNO)C(O)=O",          # N(6)-hydroxy-L-lysine (should be positive)
        "CC(C)[C@H](N(O)O)C(O)=O",         # N,N-dihydroxy-L-isoleucine (should be positive)
        "ONCC(O)=O",                     # N-hydroxyglycine (should be positive)
        "N[C@@H](CCCC)C(O)=O",            # Amino acid lacking N-hydroxy substitution (should be negative)
        "OC(=O)NO",                      # Hydroxycarbamic acid (should be negative)
        "OC1(O)C23N(CC1)C(=[N+]([O-])C(C3N=C(N2)N)COC(=O)NO)N"  # problematic false positive from previous attempt
    ]
    for smi in test_smiles:
        result, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("--------------------------------------------------")