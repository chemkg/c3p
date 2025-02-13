"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: secondary ammonium ion
Definition: An organic cation obtained by protonation of any secondary amino compound.
A protonated secondary amine (R2NH) becomes R2NH2+, meaning that the nitrogen
should have 2 heavy (organic) substituents and 2 hydrogen atoms, and a +1 charge.
While ideally the nitrogen is sp3 and non‐aromatic, the sp3 assignment may not be reliable,
so we relax that constraint and explicitly check for the two heavy substituents and 2 hydrogens.
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule (provided as a SMILES string) belongs to the class of secondary ammonium ions.
    We are looking for at least one nitrogen atom with a formal charge +1 that has exactly two heavy (non-hydrogen)
    neighbors and exactly two hydrogens attached. Aromatic nitrogens are skipped to avoid misclassifying conjugated systems.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a secondary ammonium ion is identified, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic check: ensure the molecule contains at least one carbon (i.e. is organic)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms found – not an organic molecule"
    
    # Add explicit hydrogens so that hydrogen counts are accurate.
    mol = Chem.AddHs(mol)
    
    # Loop over all atoms to identify candidate protonated amines.
    for atom in mol.GetAtoms():
        # Look for nitrogen atoms carrying a positive formal charge.
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Skip aromatic nitrogens (likely part of conjugated systems).
            if atom.GetIsAromatic():
                continue
            
            # Count heavy (non-hydrogen) neighbor atoms.
            heavy_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1)
            # Count all attached hydrogens (explicit + implicit).
            h_count = atom.GetTotalNumHs()
            
            # For a protonated secondary amine (R2NH2+), we expect exactly two heavy neighbors and two hydrogens.
            if heavy_neighbors == 2 and h_count == 2:
                return True, (f"Found secondary ammonium ion: nitrogen with {heavy_neighbors} organic substituents "
                              f"and {h_count} hydrogens (protonated secondary amine)")
    
    return False, "No protonated secondary amine (secondary ammonium ion) found"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "C[NH2+]C",  # dimethylaminium should be a secondary ammonium ion.
        "C(C[NH3+])CCC[NH2+]CC([C@@H](CO)O)=O",  # N-(D-erythrulosyl)-cadaverine contains a candidate.
        "[C@]12([C@]3([C@@]([C@@]4(C(C[C@@H](O)CC4)=CC3)C)(CC[C@]2(C)[C@]5([C@@H]([C@]6(O[C@]5(C1)[H])CC[C@@H](C)C[NH2+]6)C)[H])[H])[H])[H]",  # solasodine
        "O=S(=O)(N1CCC[NH2+]CC1)C1=CC=CC2=C1C=CN=C2"  # fasudil
    ]
    for smi in test_smiles:
        result, reason = is_secondary_ammonium_ion(smi)
        print(smi, "->", result, "|", reason)