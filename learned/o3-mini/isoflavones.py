"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: Isoflavones – any isoflavonoid with a 3-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.
"""

from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone (i.e. an isoflavonoid with a 3-aryl-1-benzopyran-4-one skeleton)
    based on its SMILES string.
    
    The strategy is two-fold:
      1. Identify the benzopyran-4-one core (i.e. the 1-benzopyran-4-one fused bicyclic system).
      2. Verify that at least one aromatic (phenyl) ring is substituted onto the core (typically at position 3).
         This is done by locating a phenyl ring that is connected to one or more of the core atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as an isoflavone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string to create an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for the benzopyran-4-one core.
    # This pattern represents a fused bicyclic system with a lactone (the 4-one).
    core_smarts = "c1ccc2c(c1)oc(=O)c2"
    core_mol = Chem.MolFromSmarts(core_smarts)
    if core_mol is None:
        return False, "Error in constructing SMARTS for benzopyran-4-one core"
    
    # Find all matches for the core in the molecule.
    core_matches = mol.GetSubstructMatches(core_mol)
    if not core_matches:
        return False, "Molecule does not contain the 1-benzopyran-4-one (benzopyranone) core"
    
    # Combine all atom indices from all core matches into a set.
    core_atoms = set()
    for match in core_matches:
        core_atoms.update(match)
    
    # Now, define a SMARTS for a phenyl ring (aromatic benzene).
    phenyl_smarts = "c1ccccc1"
    phenyl_mol = Chem.MolFromSmarts(phenyl_smarts)
    if phenyl_mol is None:
        return False, "Error in constructing SMARTS for phenyl group"
    
    # Find all phenyl group matches.
    phenyl_matches = mol.GetSubstructMatches(phenyl_mol)
    if not phenyl_matches:
        return False, "Molecule has the benzopyran-4-one core but no phenyl (aryl) substituent found"

    # For each phenyl match, check if it is attached to the core.
    # We require that at least one atom in the phenyl ring is directly bonded to an atom of the core.
    for ph_match in phenyl_matches:
        ph_set = set(ph_match)
        # Iterate over atoms in this phenyl match.
        for idx in ph_set:
            atom = mol.GetAtomWithIdx(idx)
            # Check neighbors of this atom.
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # If a neighbor is in the core atoms but is NOT itself part of the phenyl ring,
                # then we consider this phenyl group attached to the core.
                if nbr_idx in core_atoms and nbr_idx not in ph_set:
                    return True, "Molecule contains the benzopyran-4-one core with an attached aryl substituent"
    
    # If no appropriate connection is found then the molecule lacks the required 3-aryl substitution.
    return False, "Molecule contains the benzopyran-4-one core but lacks an attached aryl substituent (3-aryl group)"

# Example usage (can be used for testing):
if __name__ == "__main__":
    # Test with daidzein (an isoflavone). Its SMILES:
    test_smiles = "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O"
    result, reason = is_isoflavones(test_smiles)
    print("Result:", result)
    print("Reason:", reason)