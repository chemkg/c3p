"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
#!/usr/bin/env python3
"""
Classifies: N-acyl-L-alpha-amino acid
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An N-acyl-L-alpha-amino acid should have an L-alpha amino acid core
(i.e. a single chiral alpha-carbon attached to a carboxyl group) and at least one 
nitrogen (which may be the alpha amine or a side-chain amine, e.g. in lysine) that is 
acylated (i.e. bound to a carbonyl group via an amide bond)—but not part of a peptide chain.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    
    The algorithm works in two steps:
      1. It finds an L-alpha-amino acid core by searching for a unique chiral carbon
         (using [C@H] or [C@@H]) bonded to a carboxyl group (C(=O)O or C(=O)[O-]).
         If there is not exactly one such center, we assume the molecule is either not a 
         single amino acid or is a peptide (and thus is not classified).
      2. It then checks for the existence of an acylated nitrogen attached (directly or via a short chain)
         to that amino acid. An acylated nitrogen is defined here as a nitrogen that is bonded to a 
         carbon (other than the carboxyl carbon of the amino acid) which in turn is double‐bonded to an oxygen.
         We allow a topological distance of up to 5 bonds between the alpha‐carbon and the acylated nitrogen so that
         we can capture both direct (alpha‑N‑acyl) and side‐chain (e.g. lysine epsilon‑N‑acyl) modifications.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an N-acyl-L-alpha-amino acid, False otherwise.
        str: Explanation of the reasoning.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Step 1: Identify the L-alpha-amino acid core.
    # We build SMARTS patterns for a carboxylated chiral alpha-carbon.
    # We look for protonated (C(=O)O) and deprotonated (C(=O)[O-]) versions.
    patterns = [
        Chem.MolFromSmarts("[C@H](C(=O)O)"), 
        Chem.MolFromSmarts("[C@H](C(=O)[O-])"),
        Chem.MolFromSmarts("[C@@H](C(=O)O)"), 
        Chem.MolFromSmarts("[C@@H](C(=O)[O-])")
    ]
    alpha_indices = set()
    for pat in patterns:
        matches = mol.GetSubstructMatches(pat)
        # The first atom in our SMARTS is the chiral alpha-carbon
        for match in matches:
            alpha_indices.add(match[0])
    
    if len(alpha_indices) != 1:
        return False, ("Found %d potential L-alpha-amino acid cores; expected exactly 1. "
                       "The molecule may be a peptide or not a simple amino acid." % len(alpha_indices))
    
    alpha_idx = list(alpha_indices)[0]
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    
    # Identify the carboxyl carbon attached to the alpha carbon.
    # (Assumes that one of the neighbors of alpha is a carbon part of a carboxyl group.)
    acid_carbon_idx = None
    for nbr in alpha_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # carbon
            # Check if neighbor has at least one oxygen connected via a double bond.
            for nn in nbr.GetNeighbors():
                if nn.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        acid_carbon_idx = nbr.GetIdx()
                        break
            if acid_carbon_idx is not None:
                break
    if acid_carbon_idx is None:
        return False, "No carboxyl group found attached to the L-alpha-carbon."

    # --- Step 2: Search for an N-acyl substituent.
    # Define a helper function to decide whether a given nitrogen is acylated,
    # meaning that it is bonded to at least one carbon (other than the identified carboxyl carbon)
    # that bears a double-bonded oxygen.
    def is_acylated_nitrogen(n_atom):
        for nbr in n_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != acid_carbon_idx:
                # Check if this carbon has a double-bonded oxygen neighbor.
                for nn in nbr.GetNeighbors():
                    if nn.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                        if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            return True
        return False

    # Set a maximum topological distance from the alpha-carbon to the acylated nitrogen.
    max_distance = 5
    found_acyl = False
    # Iterate over all nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            path = rdmolops.GetShortestPath(mol, alpha_idx, atom.GetIdx())
            if 0 < len(path) - 1 <= max_distance:
                if is_acylated_nitrogen(atom):
                    found_acyl = True
                    break

    if found_acyl:
        return True, ("Molecule has an L-alpha-amino acid core (alpha carbon index %d with carboxyl group) "
                      "and an acylated nitrogen (within %d bonds) attached to it." % (alpha_idx, max_distance))
    else:
        return False, "No N-acyl substituent found attached to the L-alpha-amino acid moiety."

# Example testing (uncomment to run local tests)
# test_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"  # N-acetyl-L-aspartic acid (expected: True)
# result, reason = is_N_acyl_L_alpha_amino_acid(test_smiles)
# print(result, reason)