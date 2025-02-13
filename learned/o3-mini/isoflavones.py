"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: Isoflavones – any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.
"""

from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone 
    (i.e. an isoflavonoid with a 3-aryl-1-benzopyran-4-one skeleton)
    based on its SMILES string.

    Our two-part strategy is:
      1. Identify a benzopyran-4-one (chromen-4-one) core.
         Here we use a SMARTS that is a bit more permissive than the previous attempt.
         The pattern "c1coc2c(c1=O)cccc2" should capture the bicyclic core even if substituted.
      2. Identify an attached phenyl (aromatic benzene) that is linked to the core via a non-ring bond.
         The rationale is that an isoflavone has a 3-aryl substituent attached by a single bond 
         (as opposed to being fused into the core). 

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an isoflavone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more permissive SMARTS for the benzopyran-4-one core.
    # This pattern is intended to capture a 1-benzopyran-4-one (chromen-4-one) scaffold
    # even if extra substituents are present on the aromatic rings.
    core_smarts = "c1coc2c(c1=O)cccc2"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error constructing SMARTS for benzopyran-4-one core"
    
    # Find matches of the core in the molecule.
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Molecule does not contain the benzopyran-4-one (chromen-4-one) core"
    
    # Collect all atom indices that are in some core match.
    core_atoms = set()
    for match in core_matches:
        core_atoms.update(match)
    
    # Define a SMARTS for a phenyl (benzene) ring.
    phenyl_smarts = "c1ccccc1"
    phenyl_query = Chem.MolFromSmarts(phenyl_smarts)
    if phenyl_query is None:
        return False, "Error constructing SMARTS for phenyl group"
    
    # Find all phenyl substructure matches.
    phenyl_matches = mol.GetSubstructMatches(phenyl_query)
    if not phenyl_matches:
        return False, "Molecule has the core but no phenyl (aryl) substituent found"
    
    # Check if any phenyl ring is attached to the benzopyran-4-one core by a single (non-fused) bond.
    # We iterate over each phenyl match and for each atom in the phenyl, we look for a neighbor that is in the core.
    # We also check that the bond connecting them is NOT in a ring.
    for phenyl in phenyl_matches:
        phenyl_set = set(phenyl)
        for idx in phenyl_set:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Look for a neighbor that is in the core (but not part of the phenyl ring itself)
                if nbr_idx in core_atoms and nbr_idx not in phenyl_set:
                    bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                    # Accept only if the attaching bond is not fused (i.e. not in any ring)
                    if bond is not None and not bond.IsInRing():
                        return True, "Molecule contains the benzopyran-4-one core with an appropriately attached (non-fused) aryl substituent"
    
    return False, "Molecule contains the benzopyran-4-one core but lacks a properly attached (non-fused) 3-aryl substituent"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with daidzein (an isoflavone) SMILES
    test_smiles = "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O"
    result, reason = is_isoflavones(test_smiles)
    print("Result:", result)
    print("Reason:", reason)