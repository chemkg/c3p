"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: any organic heteromonocyclic compound with a structure based on a dihydropyrrole (pyrroline)
Definition: A pyrroline derivative (dihydropyrrole) has a five‐membered ring containing exactly one nitrogen and four carbons, 
with exactly one double bond among the bonds that connect the ring atoms.
This heuristic forces Kekulization (to reveal explicit bond orders) and uses RDKit’s GetSymmSSSR() function 
to capture candidate rings even in fused systems.
Note: This approach is heuristic and further refinements are possible.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline derivative (i.e. based on a dihydropyrrole ring)
    by checking for a five‐membered ring that contains exactly one nitrogen and four carbons,
    with exactly one double bond (based on kekulized bond orders) within that ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a suitable dihydropyrrole ring is found, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Require that molecule is organic (has at least one carbon atom).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule lacks carbon atoms; not organic."
        
    # Force Kekulization to get explicit bond orders.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        return False, "Kekulization failed; ambiguous aromaticity."
        
    # Get rings using the symmetric SSSR, which can capture rings even in fused systems.
    sssr = Chem.GetSymmSSSR(mol)
    
    # Iterate over each ring candidate.
    for ring in sssr:
        atom_indices = list(ring)  # convert the ring object (tuple-like) to a list of atom indices
        if len(atom_indices) != 5:
            continue  # consider only five-membered rings
        
        # Count number of nitrogen and carbon atoms in the ring.
        n_count = 0
        c_count = 0
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                n_count += 1
            elif atom.GetAtomicNum() == 6:
                c_count += 1
        # Require exactly one nitrogen and four carbons.
        if n_count != 1 or c_count != 4:
            continue
        
        # Get the bond indices that belong to this ring.
        bond_indices = ring.BondIndices()
        double_bond_count = 0
        for bidx in bond_indices:
            bond = mol.GetBondWithIdx(bidx)
            # Count explicit double bonds from the kekulized molecule.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1

        # We require exactly one double bond in the ring.
        if double_bond_count == 1:
            return True, "Found a five‐membered dihydropyrrole ring (1 nitrogen, 4 carbons, 1 double bond)."
    
    # If no candidate ring satisfies the criteria.
    return False, "No valid five‐membered dihydropyrrole ring (1 nitrogen, 4 carbons, and 1 double bond) found."

# Example usage for internal testing.
if __name__ == "__main__":
    test_smiles = [
        "O=C/1NCC(\\C1=C(\\O)/C=C/C(=C/C#C/C=C/C)/C)=O",  # Ravynic acid (true pyrroline derivative)
        "S=C1NCCC1",                                    # Pyrrolidine-2-thione (false: lacks the required double bond)
        "C1(N=CCC1)(C)C",                               # 5,5-dimethyl-1-pyrroline (true)
        "OC1=C2[N+](=CC=C1)[C@@H](/C=C(/C=C/C3CC3)\\C)[C@@H]([C@]2(O)C)O",  # Cyclizidine F (true pyrroline derivative)
        "O=C1N(C(O)(CC)C(=C1C(=O)C(CCCCCCC)C)O)C",      # Penicillenol D (true)
        "O=C1N[C@](C(=O)N)(CC(C)C)[C@](C1)(O)C",         # Monascustin (false by our definition)
    ]
    for smi in test_smiles:
        result, reason = is_pyrroline(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*40}")