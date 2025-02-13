"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: Primary Amine
Definition: A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.
The program looks for an –NH2 group (i.e. primary amine) that is attached to at least one carbon and is not part of an amide.
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine has an –NH2 group where nitrogen is bound to one hydrocarbyl substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a primary amine, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a nitrogen with exactly two attached hydrogens (i.e. potentially an NH2 group).
    # This pattern will match both aliphatic and aromatic primary amines.
    primary_amine_smarts = "[NX3;H2]"
    primary_amine_pattern = Chem.MolFromSmarts(primary_amine_smarts)
    if primary_amine_pattern is None:
        return False, "Error in generating SMARTS pattern"
    
    # Find all matches of the primary amine pattern in the molecule.
    matches = mol.GetSubstructMatches(primary_amine_pattern)
    if not matches:
        return False, "No primary amine (-NH2) substructure found"
    
    # Check each matched nitrogen for required features:
    # It must be attached to at least one carbon atom (i.e. the hydrocarbyl group),
    # and it must not be part of an amide (i.e. attached to a carbonyl carbon, C(=O)-).
    for match in matches:
        # match is a tuple of atom indices that match the pattern.
        # Our pattern only covers one atom (the nitrogen) so take the first index.
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Flag to check if this nitrogen is attached to at least one carbon atom.
        attached_to_carbon = False
        amide_flag = False
        
        for neighbor in n_atom.GetNeighbors():
            # Check if the neighbor atom is carbon (atomic number 6) to represent a hydrocarbyl group.
            if neighbor.GetAtomicNum() == 6:
                attached_to_carbon = True
                # Further check if this carbon is part of a carbonyl (C=O) bond,
                # which would indicate an amide group.
                for bond in neighbor.GetBonds():
                    # If the bond is to an oxygen and is a double bond, then it is likely a carbonyl.
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                        amide_flag = True
                        break
            # If we already detected an amide connection for this candidate, break early.
            if amide_flag:
                break

        if amide_flag:
            # This matched amine group appears to be part of an amide.
            continue
        
        if not attached_to_carbon:
            # If not bound to any carbon, it is not derived from ammonia by replacing a hydrogen with a hydrocarbyl group.
            continue
        
        # Found a valid primary amine substructure.
        return True, "Molecule contains a primary amine (R–NH2) group."
    
    return False, "No valid primary amine group (bound to a hydrocarbyl moiety and not part of an amide) found."