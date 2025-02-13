"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: Carbonate Ester
Definition: Any carbonate that is carbonic acid in which at least one (or both) of the –OH groups is (are)
replaced by organyl groups. This includes diesters (e.g. dimethyl carbonate) and monoesters (e.g. monoethyl carbonate)
as well as cyclic carbonate derivatives.
"""

from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester is defined as a derivative of carbonic acid (O=C(OH)2)
    in which at least one hydroxyl (-OH) is replaced by an organyl group.
    That is, we look for a carbon atom that is connected via a double bond to an oxygen,
    and via two single bonds to two oxygens—with no direct carbon attachments!
    Finally, at least one of the oxygen substituents (apart from the carbonyl oxygen)
    must be bound to a non-hydrogen atom (indicating an ester substitution).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we can check for -OH groups and substitutions
    mol = Chem.AddHs(mol)
    
    # We now search for a carbon atom with the carbonate connectivity:
    # * The carbon (C_carbonate) must be bonded to a double-bonded oxygen (the carbonyl oxygen)
    # * C_carbonate must also be bonded via single bonds to two oxygens.
    # * C_carbonate must have no additional (direct) substituents (i.e. no direct carbon attachment).
    # This distinguishes the carbonate motif (O–C(=O)–O) from a normal carboxylic ester.
    for atom in mol.GetAtoms():
        # Look at carbon atoms only.
        if atom.GetAtomicNum() != 6:
            continue
        neighbors = atom.GetNeighbors()
        
        # For a valid carbonate moiety the carbon should have exactly 3 neighbors
        if len(neighbors) != 3:
            continue
        
        double_bonded_oxygens = []
        single_bonded_oxygens = []
        is_valid_connectivity = True
        
        # Examine bonds to the carbon atom
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            # We are only interested in oxygen neighbors.
            if nbr.GetAtomicNum() != 8:
                is_valid_connectivity = False
                break
            # Classify the bond type.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                double_bonded_oxygens.append(nbr)
            elif bond.GetBondType() == Chem.BondType.SINGLE:
                single_bonded_oxygens.append(nbr)
            else:
                # If bond type is not single or double, ignore this candidate.
                is_valid_connectivity = False
                break
        if not is_valid_connectivity:
            continue
        
        # Check for exactly one double-bonded oxygen and exactly two single-bonded oxygens.
        if len(double_bonded_oxygens) != 1 or len(single_bonded_oxygens) != 2:
            continue
        
        # At this point we have a "carbonate skeleton" (O–C(=O)–O).
        # Now check the substitution on the two single-bonded oxygens.
        # In carbonic acid, both oxygens are -OH (each oxygen is attached to a hydrogen).
        # In a carbonate ester, at least one of these oxygens is bound to a non-hydrogen (an organyl group).
        substituted_found = False
        for oxy in single_bonded_oxygens:
            # List the neighbors of this oxygen (excluding the candidate central carbon).
            other_neighbors = [nbr for nbr in oxy.GetNeighbors() if nbr.GetIdx() != atom.GetIdx()]
            # In a typical hydroxyl group the only other neighbor is a hydrogen.
            # If any of these neighbors is not hydrogen (atomic number != 1), we say that oxygen is substituted.
            if any(nbr.GetAtomicNum() != 1 for nbr in other_neighbors):
                substituted_found = True
                # We do not require both to be substituted:
                # a monoester counts as a carbonate ester.
                break

        if substituted_found:
            return True, "Found carbonate moiety (O–C(=O)–O) with at least one organyl substitution."
        else:
            # If we find the carbonate skeleton but both oxygens are -OH then it is just carbonic acid.
            return False, "Found a carbonic acid moiety (O=C(OH)2) without any ester substitution."
    
    return False, "No carbonate ester motif (O–C(=O)–O) detected in the molecule."

# For quick testing (uncomment to run)
# test_smiles = [
#     "COC(=O)OC",  # dimethyl carbonate
#     "CCOC(O)=O",  # monoethyl carbonate
#     "O1CCOC1=O", # ethylene carbonate (cyclic)
#     "O=C(O)O"    # carbonic acid (should be False)
# ]
# for smi in test_smiles:
#     result, reason = is_carbonate_ester(smi)
#     print(f"SMILES: {smi} --> {result}. Reason: {reason}")