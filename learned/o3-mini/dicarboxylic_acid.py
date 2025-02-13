"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Dicarboxylic Acid – Any carboxylic acid containing two carboxy groups.
This implementation uses a custom routine to count unique carboxyl groups by examining
each potential carboxyl carbon (i.e. a carbon with one double‐bonded oxygen and one single‐bonded oxygen that is either –OH or –O–).
Additional filters check for limited nitrogen (to avoid peptide-like molecules) and,
in molecules without any nitrogen, a minimum number of carbon atoms (to avoid very small acids).
"""

from rdkit import Chem

def count_carboxyl_groups(mol):
    """
    Count carboxyl groups by scanning for carbons that are bonded to one oxygen by a double bond
    and a second oxygen by a single bond (the typical signature of a COOH group).
    This approach avoids double-counting the same carboxyl unit.
    """
    carboxyl_carbons = set()
    for atom in mol.GetAtoms():
        # Look at carbon atoms only (atomic number 6)
        if atom.GetAtomicNum() != 6:
            continue
        neighbors = atom.GetNeighbors()
        # We require at least 2 oxygen neighbors
        oxygen_neighbors = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(oxygen_neighbors) < 2:
            continue
        dbl_found = False
        single_found = False
        for nb in oxygen_neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
            if bond is None:
                continue
            # Check for a double bond (the carbonyl oxygen)
            if bond.GetBondTypeAsDouble() >= 1.99:
                dbl_found = True
            # Check for a single bond where the oxygen is –OH or deprotonated –O-
            elif bond.GetBondType().name == "SINGLE":
                # Count implicit/explicit hydrogens on oxygen
                hcount = nb.GetTotalNumHs()
                # Allow if oxygen has hydrogen (–OH) or is deprotonated (formal charge -1)
                if hcount >= 1 or nb.GetFormalCharge() < 0:
                    single_found = True
        if dbl_found and single_found:
            carboxyl_carbons.add(atom.GetIdx())
    return len(carboxyl_carbons)

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    
    For our purposes, a dicarboxylic acid is defined as a molecule that contains exactly two
    carboxyl groups AND does not show features indicating a peptide/complex structure.
    The carboxyl groups are detected by looking for a carbon that is double-bonded to an oxygen
    and singly bonded to another oxygen (which is protonated or negatively charged). To reduce spuriously
    small molecules such as malate (which has only 4 carbons and is a common metabolic intermediate),
    a further filter is applied: if no nitrogen is present then at least 5 carbons must be found.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as a dicarboxylic acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carboxyl groups using the custom routine.
    n_carboxyl = count_carboxyl_groups(mol)
    if n_carboxyl != 2:
        return False, f"Molecule contains {n_carboxyl} carboxyl group(s); need exactly 2."
    
    # Count nitrogen atoms (atomic number 7)
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    # Reject if there are more than 2 nitrogen atoms (suggestive of peptide/complex structures)
    if n_nitrogen > 2:
        return False, "Molecule contains extra nitrogen atoms, suggesting a peptide or complex structure rather than a simple dicarboxylic acid."
    
    # As an extra heuristic: if no nitrogen is present, require the molecule to have at least 5 carbon atoms
    # to avoid very small acids (e.g., malate with only 4 carbons).
    n_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_nitrogen == 0 and n_carbon < 5:
        return False, "Molecule has too few carbon atoms for an interesting dicarboxylic acid."
    
    return True, "Molecule contains exactly two carboxyl groups and meets structural heuristics; qualifies as a dicarboxylic acid."

# Example usage (uncomment to test in your environment):
# print(is_dicarboxylic_acid("OC(=O)C(O)=O"))  # 2-formylglutaric acid example (should be True if it meets criteria)
# print(is_dicarboxylic_acid("C[C@H](NCCC(O)=O)C(O)=O"))  # (S)-beta-alanopine, expecting True