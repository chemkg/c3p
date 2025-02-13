"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: Aldoxime (Oximes of aldehydes RCH=NOH)

An aldoxime is formed from an aldehyde by conversion of the carbonyl group (RCH=O) to an oxime group (RCH=NOH).
Thus, the structure must feature a carbon–nitrogen double bond (C=N) with the nitrogen bonded to a hydroxyl (-OH)
group, and importantly the carbon must have exactly one hydrogen.

This script converts a SMILES string to a molecule (with explicit hydrogens added) and then checks for the following:
  1. A nitrogen atom that is double-bonded to a carbon that has exactly one hydrogen.
  2. That same nitrogen is single-bonded to an oxygen that carries a hydrogen (i.e. an –OH group).
If such a set of atoms is found, the molecule is classified as an aldoxime.
"""

from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime (RCH=NOH) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an aldoxime, False otherwise.
        str: Reason for the classification.
    """
    # Attempt to create a molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are correct.
    mol = Chem.AddHs(mol)
    
    # Iterate over all nitrogen atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Only interested in nitrogen
        
        # Initialize variables to check for required bonds.
        carbon_candidate = None
        oxygen_candidate = None
        
        # Iterate over neighbors of the nitrogen to find a carbon (double bond) and an oxygen (single bond with -OH).
        for neighbor in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            # Check for oxygen neighbor (for -OH group)
            if neighbor.GetAtomicNum() == 8:
                # For the oxime, the oxygen should be connected via a single bond and have exactly one hydrogen.
                if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetTotalNumHs() == 1:
                    oxygen_candidate = neighbor
            # Check for carbon neighbor (the aldehyde carbon) connected via a double bond.
            elif neighbor.GetAtomicNum() == 6:
                if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetTotalNumHs() == 1:
                    carbon_candidate = neighbor
        
        # If both an appropriate carbon and oxygen candidate were found, we identify the aldoxime group.
        if carbon_candidate is not None and oxygen_candidate is not None:
            return True, "Contains an aldoxime group (RCH=NOH)"
    
    # If no valid aldoxime group is found.
    return False, "Aldoxime group not found in the molecule"

# Uncomment the code below for simple testing examples:
# test_smiles = [
#     "[H]\\C(C(C)C)=N/O",                     # (E)-2-methylpropanal oxime
#     "C([C@@H](/C(=N/O)/[H])C)C",             # (1E,2S)-2-methylbutanal oxime
#     "COc1cc(\\C=N/O)nc(c1)-c1ccccn1",          # (Z)-4-methoxy-2,2-bipyridine-6-carbaldehyde oxime
#     "COc1cc(\\C=N\\O)nc(-c2ccccn2)c1OC",       # caerulomycin C
#     "C(C=NO)=NO",                            # glyoxime
#     "[H]C(C)=NO",                           # acetaldehyde oxime
#     "CN1C=CCC=C1\\C=N\\O"                     # ProPAM
# ]
# 
# for smi in test_smiles:
#     result, reason = is_aldoxime(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")