"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: nitrile compounds (RC≡N)
A nitrile is defined as a compound featuring a carbon–nitrogen triple bond 
where the carbon is substituted (i.e. not simply HC≡N). 
This version adds explicit hydrogens and then for each triple bond between 
a carbon and a nitrogen, it checks that the nitrile carbon has exactly exactly two neighbors:
one must be the nitrile nitrogen and the other must be a non‐hydrogen substituent.
Metal substituents are accepted.
"""

from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a substituted nitrile (RC≡N) based on its SMILES string.
    
    The method adds explicit hydrogens so that the bond environment is unambiguous.
    It then iterates over all triple bonds between a carbon and a nitrogen.
    For the nitrile carbon the algorithm requires that it is connected to exactly two atoms:
    (a) the nitrile nitrogen and (b) one substituent; if that substituent is a hydrogen,
    then the motif is simply HCN and not RC≡N.
    Metal substituents are now accepted.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains at least one valid substituted nitrile group.
        str: An explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the degree count is unambiguous.
    mol = Chem.AddHs(mol)
    
    # Iterate over all bonds in the molecule.
    for bond in mol.GetBonds():
        # Check for a triple bond.
        if bond.GetBondType() != Chem.BondType.TRIPLE:
            continue
        
        # Identify the two atoms in the bond.
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        
        # One atom must be carbon and the other nitrogen.
        if (a1.GetSymbol() == "C" and a2.GetSymbol() == "N"):
            nitrile_c = a1
            nitrile_n = a2
        elif (a2.GetSymbol() == "C" and a1.GetSymbol() == "N"):
            nitrile_c = a2
            nitrile_n = a1
        else:
            continue
        
        # After explicitly adding hydrogens, a C in a proper RC≡N group (substituted nitrile)
        # should have exactly 2 neighbors: the nitrile N and one substituent.
        neighbors = nitrile_c.GetNeighbors()
        if len(neighbors) != 2:
            continue
        
        # Identify the substituent that is not the nitrile nitrogen.
        substituent = None
        for nb in neighbors:
            if nb.GetIdx() == nitrile_n.GetIdx():
                continue
            substituent = nb
            break
        
        if substituent is None:
            continue
        
        # If the substituent is a hydrogen, then it is simply HC≡N (not a substituted nitrile).
        if substituent.GetAtomicNum() == 1:
            continue
        
        # If we reached here, we found a nitrile carbon (C with a triple bond to N) that is 
        # connected only to the nitrile nitrogen and to a substituent that is not hydrogen. 
        # This qualifies as RC≡N.
        return True, "Molecule contains a substituted nitrile group (RC≡N)"
    
    # If no valid nitrile groups were found, report the reason.
    return False, "No substituted nitrile group (RC≡N) found; only HCN‐like or invalid nitrile(s) detected"

# Example testing if run as script:
if __name__ == "__main__":
    test_smiles_list = [
        "CC#N",                    # Acetonitrile: valid nitrile (substituted with a CH3)
        "C#N",                     # Hydrogen cyanide: should be rejected
        "N#C[Fe]C#N",              # Iron dicyanide: valid by our new rules (metal substituent accepted)
    ]
    for smi in test_smiles_list:
        result, reason = is_nitrile(smi)
        print(f"{smi}: {result} -- {reason}")