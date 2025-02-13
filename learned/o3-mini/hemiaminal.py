"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: Hemiaminal
Definition: “Any organic amino compound that has an amino group and a hydroxy group attached 
to the same carbon atom. Hemiaminals are intermediates in the formation of imines by addition of 
an amine to an aldehyde or ketone; those derived from primary amines are particularly unstable.”
 
This version adds explicit hydrogens and, for every sp3 carbon with four bonds, it checks that:
    - One neighbor is an oxygen with exactly two neighbors (the candidate carbon itself and one hydrogen)
      to indicate a genuine -OH group.
    - One neighbor is a nitrogen that (a) is connected by a single bond and (b) does not show 
      any sign (i.e. double bond) of being part of an amide.
If both conditions are met (for at least one candidate carbon) the molecule is declared to contain 
a hemiaminal functional group.
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has a tetrahedral (mostly sp3) carbon that bears both a hydroxyl group (-OH)
    and an amino group (-NH2, -NHR, or -NR2) attached directly to that same carbon.
    
    The approach is:
      1. Parse the molecule and add explicit hydrogens.
      2. Iterate over sp3 carbon atoms that have 4 bonds.
      3. For each such carbon, check that one neighbor is a hydroxyl oxygen defined
         by having exactly 2 neighbors (the candidate carbon and one hydrogen).
      4. Also check that one neighbor is a nitrogen that does not have any double bond
         to an oxygen (which might indicate an amide).
      5. If a carbon matches both criteria, report success and return the indices.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a hemiaminal is detected, False otherwise.
        str: Reason/message with matching atom indices or explanation.
    """
    # Try to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string; could not parse molecule"
        
    # Add explicit hydrogens to help detect -OH groups.
    mol = Chem.AddHs(mol)
    
    hemiaminal_matches = []  # we will collect tuples (carbon_idx, oxygen_idx, nitrogen_idx)
    
    # Helper function: Is the oxygen a genuine -OH group?
    def is_hydroxyl(oxygen, parent_idx):
        # Check that the oxygen is connected to the candidate carbon (parent_idx)
        # and has exactly one hydrogen and no extra heavy atom
        nbrs = oxygen.GetNeighbors()
        # After AddHs, oxygen in a -OH should have exactly 2 neighbors: the carbon and one hydrogen.
        if len(nbrs) != 2:
            return False
        has_h = False
        valid = False
        for nbr in nbrs:
            # One neighbor must be a hydrogen.
            if nbr.GetAtomicNum() == 1:
                has_h = True
            else:
                # The only heavy atom neighbor must be our candidate carbon.
                if nbr.GetIdx() != parent_idx:
                    return False
                else:
                    valid = True
        return has_h and valid

    # Helper function: Is the nitrogen a free (non-amide) amino group?
    def is_valid_amino(nitrogen):
        # Ensure no bond from nitrogen is a double bond to an oxygen.
        for nbr in nitrogen.GetNeighbors():
            # Get the connecting bond.
            bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), nbr.GetIdx())
            if bond is None:
                continue
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                return False
        return True

    # Iterate over all atoms to search for sp3 carbons with 4 substituents.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # Only consider carbons.
        # Must be sp3 hybridized.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # After adding all explicit hydrogens, the total number of bonds should be exactly 4.
        if atom.GetDegree() != 4:
            continue
        
        current_idx = atom.GetIdx()
        o_matches = []  # list of oxygen neighbor indices that pass the -OH check.
        n_matches = []  # list of nitrogen neighbor indices that pass the amino check.
        
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(current_idx, nbr.GetIdx())
            # Only consider single bonds.
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            atomic_num = nbr.GetAtomicNum()
            if atomic_num == 8:
                # Check that the oxygen is in a genuine -OH group.
                if is_hydroxyl(nbr, current_idx):
                    o_matches.append(nbr.GetIdx())
            elif atomic_num == 7:
                # Check for amide character: no double bonded oxygen.
                if is_valid_amino(nbr):
                    n_matches.append(nbr.GetIdx())
        
        # Record a match if at least one valid -OH and one valid -NHx are present.
        if o_matches and n_matches:
            for o_idx in o_matches:
                for n_idx in n_matches:
                    hemiaminal_matches.append((current_idx, o_idx, n_idx))
    
    if hemiaminal_matches:
        return True, f"Found hemiaminal substructure in atoms with indices: {hemiaminal_matches}"
    else:
        return False, "No carbon with both hydroxyl and amino substituents found; not a hemiaminal"


# Example usage
if __name__ == "__main__":
    # Test a few SMILES strings.
    test_smiles_list = [
        "OC(N)CC",  # 2-Aminopropanol: expected to be hemiaminal.
        "C1[C@@]2(N3CC[C@@]42[C@]5(N(C=6C4=CC=CC6)C(C[C@]7([C@@]5([C@@]1(C(=CCO7)C3)[H])[H])[H])=O)[H])O"  # pseudostrychnine, one of the positives.
    ]
    for s in test_smiles_list:
        result, reason = is_hemiaminal(s)
        print(result, reason)