"""
Classifies: CHEBI:46895 lipopeptide
"""
#!/usr/bin/env python3
"""
Classifies: Lipopeptide
Definition: A compound consisting of a peptide with an attached lipid.
The algorithm first checks for the presence of at least one amide bond as a marker for the peptide portion.
Then it computes the longest contiguous chain of non‐aromatic carbon atoms (using a DFS).
If that chain has at least eight carbons, it is taken as evidence for a lipid tail.
"""

from rdkit import Chem

def longest_chain_from(atom_idx, mol, visited):
    """Recursively compute longest chain length starting from atom_idx.
    Only considers neighboring atoms that are carbon (atomic number 6)
    and not aromatic. The bond type is allowed if it is SINGLE or DOUBLE.
    """
    visited.add(atom_idx)
    current_length = 1
    max_extension = 0
    atom = mol.GetAtomWithIdx(atom_idx)
    for bond in atom.GetBonds():
        # Get the neighboring atom's index
        nbr = bond.GetOtherAtom(atom).GetIdx()
        nbr_atom = mol.GetAtomWithIdx(nbr)
        # Only consider neighbor if it is a carbon (atomic num 6) and non‐aromatic.
        if nbr not in visited and nbr_atom.GetAtomicNum() == 6 and not nbr_atom.GetIsAromatic():
            # Allow single or double bonds (lipid chains sometimes contain unsaturation)
            if bond.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                # Use a copy of visited to allow different DFS paths
                extension = longest_chain_from(nbr, mol, visited.copy())
                if extension > max_extension:
                    max_extension = extension
    return current_length + max_extension

def get_longest_aliphatic_chain(mol):
    """Loops over all carbon atoms in the molecule (non‐aromatic)
    and returns the maximum chain length found.
    """
    max_chain = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            chain_length = longest_chain_from(atom.GetIdx(), mol, set())
            if chain_length > max_chain:
                max_chain = chain_length
    return max_chain

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide should contain a peptide portion (indicated by amide bonds)
    and an attached lipid chain (modeled as a contiguous chain of at least 8 carbons).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a lipopeptide, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES to molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of amide bonds (the peptide bond is typically C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = mol.GetSubstructMatches(amide_pattern)
    if not peptide_matches:
        return False, "No amide (peptide) bonds found"

    # Compute the longest contiguous chain of non-aromatic carbons (heuristic for lipid chain)
    longest_chain = get_longest_aliphatic_chain(mol)
    if longest_chain < 8:
        return False, f"Longest contiguous non-aromatic carbon chain is {longest_chain}, which is insufficient (need >= 8)"
    
    return True, f"Molecule contains peptide bonds and a carbon chain of length {longest_chain}"

# Example usage (for testing purposes)
if __name__ == '__main__':
    # Here we test using one of the lipopeptide examples.
    # For example, surfactin A is known to be a lipopeptide.
    surfactin_A = "[H][C@@]1(CCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1"
    result, reason = is_lipopeptide(surfactin_A)
    print("Surfactin A classification:", result, reason)