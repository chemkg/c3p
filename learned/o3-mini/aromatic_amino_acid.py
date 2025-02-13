"""
Classifies: CHEBI:33856 aromatic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: Aromatic Amino Acid
Definition: An amino acid whose structure includes an aromatic ring.
Improved criteria:
  - Must be a small molecule (heavy atom count <= 30) to avoid peptides.
  - Must have an amino acid backbone fragment of the form:
       [NX3;!$(NC(=O))][C;H](C(=O)[O;H1,-])
    (i.e. an amino group – which might be N‐substituted – directly attached to an α–carbon that bears a carboxyl or carboxylate).
  - The α–carbon must have at least one side‐chain neighbor (other than the backbone N and carboxyl carbon)
    that leads to an aromatic ring by exactly one intervening atom (i.e. a CH2 linker).
The idea is that the “proper” aromatic amino acids (e.g. phenylalanine, tyrosine, tryptophan, histidine)
have a side chain where the aromatic system is attached directly via one extra carbon.
Molecules that have a more remote aromatic ring (or additional rings) are then filtered out.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid.
    
    For classification, the molecule must:
      - Parse correctly.
      - Be small (heavy atom count <= 30).
      - Contain an amino acid backbone fragment of the form:
           [NX3;!$(NC(=O))][C;H](C(=O)[O;H1,-])
        which captures both free and N–substituted amino groups.
      - The α–carbon (the middle atom of the backbone match) must have at least one side-chain 
        neighbor (i.e. not the backbone N nor the carboxyl carbon) that is connected (by one bond)
        to an aromatic atom (i.e. typical of the common CH2–aromatic linkage).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for an aromatic amino acid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add hydrogens to have proper count for the backbone (if needed)
    mol = Chem.AddHs(mol)
    
    # Filter out molecules that are too large to be a single amino acid.
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms > 30:
        return False, f"Molecule is too large (has {heavy_atoms} heavy atoms) to be a single amino acid"
    
    # SMARTS for a generic amino acid backbone:
    #   [NX3;!$(NC(=O))]   : A trivalent nitrogen that is NOT linked to a carbonyl (to avoid amides)
    #   [C;H]              : An α–carbon with at least one hydrogen (typical for non–glycine residues)
    #   (C(=O)[O;H1,-])    : A carboxyl/carboxylate group.
    aa_smarts = "[NX3;!$(NC(=O))][C;H](C(=O)[O;H1,-])"
    aa_pattern = Chem.MolFromSmarts(aa_smarts)
    if aa_pattern is None:
        return False, "Error in parsing amino acid SMARTS pattern"
    
    # Check for the amino acid backbone substructure:
    matches = mol.GetSubstructMatches(aa_pattern)
    if not matches:
        return False, "No amino acid backbone (N–α–C–(C=O)[O]) fragment found"
    
    # For simplicity, use the first found backbone match.
    # By our SMARTS, the match indices are:
    #   match[0] : nitrogen,
    #   match[1] : α–carbon,
    #   match[2] : carboxyl carbon.
    backbone = matches[0]
    n_idx, alpha_idx, carboxyl_idx = backbone
    
    # Now ensure the aromatic ring is directly appended to the side chain.
    # The side chain atoms are neighbors of the α–carbon excluding the backbone N and carboxyl carbon.
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    sidechain_neighbors = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors() 
                           if nbr.GetIdx() not in (n_idx, carboxyl_idx)]
    
    if not sidechain_neighbors:
        return False, "No side chain attached to the α–carbon"
    
    # For each candidate side chain atom, check if one of its neighbors (other than the α–carbon)
    # is an aromatic atom (which we require to be directly attached via one bond).
    aromatic_sidechain_found = False
    for nb_idx in sidechain_neighbors:
        nb_atom = mol.GetAtomWithIdx(nb_idx)
        for nbr in nb_atom.GetNeighbors():
            # Exclude the bond back to the α–carbon.
            if nbr.GetIdx() == alpha_idx:
                continue
            # Check if this neighbor is aromatic.
            if nbr.GetIsAromatic():
                # Also, verify that the path from the α–carbon to this aromatic atom is exactly 2 bonds.
                path = rdmolops.GetShortestPath(mol, alpha_idx, nbr.GetIdx())
                if len(path) - 1 == 2:
                    aromatic_sidechain_found = True
                    break
        if aromatic_sidechain_found:
            break

    if not aromatic_sidechain_found:
        return False, "No aromatic side chain attached directly (via one intervening atom) to the α–carbon found"
    
    return True, "Molecule contains an amino acid backbone and an appropriately positioned aromatic side chain"

# Example usage (can be removed or commented out if being used as a module)
if __name__ == "__main__":
    test_examples = {
        "N-methyl-D-dopa": "CN[C@H](CC1=CC=C(O)C(O)=C1)C(O)=O",
        "D-phenylalanine": "N[C@H](Cc1ccccc1)C(O)=O",
        "N(5)-phenyl-L-glutamine (false positive expected)": "N[C@@H](CCC(=O)Nc1ccccc1)C(O)=O",
        "3-aminobenzoic acid (false negative expected)": "Nc1cccc(c1)C(O)=O"
    }
    for name, smi in test_examples.items():
        result, reason = is_aromatic_amino_acid(smi)
        print(f"{name}:\n  Result: {result}\n  Reason: {reason}\n")