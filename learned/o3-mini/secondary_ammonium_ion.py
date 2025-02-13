"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: secondary ammonium ion
Definition: An organic cation obtained by protonation of any secondary amino compound.
When a secondary amino compound (R2NH) is protonated it becomes R2NH2+.
Such nitrogen should be tetra‐coordinated (degree 4 after adding Hs) with exactly two bonds
to heavy (non‐hydrogen) atoms. Aromatic nitrogen centers are skipped.
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule (provided as a SMILES string) contains at least one
    protonated secondary amine (i.e. a secondary ammonium ion). This is defined as a
    nitrogen atom with a formal +1 charge that is tetra‐coordinated (degree 4 after H addition)
    with exactly 2 bonds to heavy (non‐hydrogen) atoms. Aromatic nitrogens are skipped.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains a protonated secondary amine, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is organic (at least one carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms found – not an organic molecule"
    
    # Add explicit hydrogens so that hydrogen counts and valence are accurate.
    mol = Chem.AddHs(mol)
    
    # Loop over all atoms to find candidate nitrogen atoms having a +1 formal charge.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Skip aromatic nitrogens (likely part of conjugated systems)
            if atom.GetIsAromatic():
                continue
            
            # Get total degree (number of bonds) from the atom including H neighbors.
            total_degree = atom.GetDegree()
            # Count heavy (non-hydrogen) neighbors.
            heavy_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1)
            
            # A protonated secondary amine (R2NH2+) should be tetra‐coordinated (degree 4)
            # with exactly 2 bonds to heavy atoms.
            if total_degree == 4 and heavy_neighbors == 2:
                return True, (f"Found protonated secondary amine (secondary ammonium ion): "
                              f"degree {total_degree} with {heavy_neighbors} heavy neighbors")
    
    return False, "No protonated secondary amine (secondary ammonium ion) found"

# Example usage (for manual testing):
if __name__ == "__main__":
    test_smiles = [
        "C[NH2+]C",  # dimethylaminium (should be positive)
        "C(C[NH3+])CCC[NH2+]CC([C@@H](CO)O)=O",  # N-(D-erythrulosyl)-cadaverine(2+); one candidate should match
        "[C@]12([C@]3([C@@]([C@@]4(C(C[C@@H](O)CC4)=CC3)C)(CC[C@]2(C)[C@]5([C@@H]([C@]6(O[C@]5(C1)[H])CC[C@@H](C)C[NH2+]6)C)[H])[H])[H])[H]",  # solasodine(1+)
        "O=S(=O)(N1CCC[NH2+]CC1)C1=CC=CC2=C1C=CN=C2",  # fasudil(1+)
        "C1CCCC[NH2+]1",  # piperidinium
        "CC[NH2+]CC"  # diethylammonium (should be secondary ammonium ion)
    ]
    for smi in test_smiles:
        result, reason = is_secondary_ammonium_ion(smi)
        print(smi, "->", result, "|", reason)