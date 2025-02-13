"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
A compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.
This program attempts to identify a nitrogen atom that is NOT quaternary (i.e. has exactly three substituents, 
where the substituents count = (number of heavy atom neighbors) + (implicit hydrogens)). 
It excludes nitrogen that is part of an “amide‐like” group – that is, a nitrogen attached (via a carbon neighbor) 
to a carbonyl group in which the carbon is double‐bonded to an oxygen (but not if that oxygen is replaced by sulfur).
We do not automatically exclude nitrogen with a positive formal charge if it still has three substituents.
"""

from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule contains an amine functional group derived from ammonia.
    
    This function scans through all nitrogen atoms in the molecule.
    For each nitrogen, it calculates the total substituents as the sum of heavy-atom neighbors
    plus the number of implicit (and explicit) hydrogens. If that sum equals 3 then the candidate may 
    be a primary (RNH2), secondary (R2NH) or tertiary (R3N) amine – even if the nitrogen bears a formal positive charge.
    
    Before labeling it as an amine, the function excludes cases in which the N is directly bonded (via a carbon)
    to a classical carbonyl group (C=O). (A C(=O)S group is allowed.)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if at least one amine group is found, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms to search for a nitrogen that fits our amine criteria.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:  # Skip non-nitrogen atoms.
            continue
        
        # Calculate the total number of substituents:
        # This equals the count of heavy-atom neighbors plus the total number of associated hydrogens.
        total_substituents = len(atom.GetNeighbors()) + atom.GetTotalNumHs()
        
        # For a primary, secondary or tertiary amine derived from ammonia, this value must be exactly 3.
        if total_substituents != 3:
            continue
        
        # Exclude nitrogen that is bonded to a carbonyl group.
        # For every neighbor of the nitrogen, if that neighbor is a carbon and is double-bonded
        # to oxygen (atomic num 8) then treat it as an amide. Here we only exclude if the oxygen is O (not S).
        is_amide = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # Carbon neighbor
                for bond in nbr.GetBonds():
                    # Check if the bond is a double bond.
                    if bond.GetBondTypeAsDouble() == 2:  # a double bond
                        other = bond.GetOtherAtom(nbr)
                        # Exclude only if the other atom is oxygen (and not sulfur)
                        if other.GetAtomicNum() == 8:
                            is_amide = True
                            break
                if is_amide:
                    break
        if is_amide:
            continue

        # If we reached here, we have found at least one nitrogen with exactly 3 substituents that is not part of a classic amide.
        return True, "Found at least one amine functional group derived from ammonia"
    
    # If no nitrogen meets our criteria then we did not detect an amine group.
    return False, "No amine functional group (derived from ammonia) identified"

# Example test run (uncomment to try)
# test_smiles = "CCCN(CCCC)CC(O)c1cc(Cl)cc2\\C(=C/c3ccc(Cl)cc3)c3cc(Cl)ccc3-c12"  # lumefantrine
# print(is_amine(test_smiles))