"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
#!/usr/bin/env python3
"""
Classifies: N-acyl-L-alpha-amino acid
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An N-acyl-L-alpha-amino acid should have an L-alpha-amino acid core
(i.e. a single chiral alpha‑carbon attached to a carboxyl group) and an acylated
nitrogen (which may be the α‑amine or a side-chain amine) that is not part of a peptide bond.
For improved performance we restrict to molecules with exactly one α‑core.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    
    The approach is as follows:
      1. Find all potential L-alpha amino acid cores by searching for a chiral carbon ([C@H] or [C@@H])
         that is directly attached to a carboxyl group (C(=O)O or C(=O)[O-]). We require exactly one such core.
         (If there are multiple cores, then the molecule probably represents a peptide rather than a single amino acid.)
      2. Identify the carboxyl carbon attached to that α‑carbon (to exclude it from being counted for acylation).
      3. Search for a nitrogen (N) atom within a topological distance (up to 5 bonds) from the α‑carbon that is “acylated.”
         A nitrogen is considered acylated if it is attached (directly or via a short chain) to a carbon (other than the carboxyl carbon)
         that bears a double‐bonded oxygen. In addition, if that acyl carbon is also attached to another nitrogen that in turn is bonded
         to an amino acid core (i.e. a peptide‐bond arrangement), then that functionality is interpreted as part of a peptide bond 
         rather than an N‑acyl substituent.
      
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
    
    # --- Step 1: Identify potential L-alpha amino acid cores.
    # We search for a chiral carbon ([C@H] or [C@@H]) directly bound to a carboxyl group.
    core_smarts = [
        "[C@H](C(=O)O)", 
        "[C@H](C(=O)[O-])",
        "[C@@H](C(=O)O)",
        "[C@@H](C(=O)[O-])"
    ]
    alpha_indices = set()
    for smarts in core_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # The first atom in the match is the chiral alpha-carbon
            alpha_indices.add(match[0])
    
    # If not exactly one L-alpha core is found, we reject.
    if len(alpha_indices) != 1:
        return False, ("Found %d potential L-alpha-amino acid cores; expected exactly 1. "
                       "The molecule may be a peptide or not a simple amino acid." % len(alpha_indices))
    
    alpha_idx = list(alpha_indices)[0]
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    
    # --- Step 2: Identify the carboxyl carbon attached to the alpha carbon.
    acid_carbon_idx = None
    for nbr in alpha_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # carbon
            # Check if this neighbor (acid carbon) bears a double-bonded oxygen.
            for nn in nbr.GetNeighbors():
                if nn.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                        acid_carbon_idx = nbr.GetIdx()
                        break
            if acid_carbon_idx is not None:
                break
    if acid_carbon_idx is None:
        return False, "No carboxyl group found attached to the L-alpha-carbon."
    
    # --- Step 3: Define helper functions.
    
    def is_peptide_bond(acyl_c_idx, candidate_n_idx):
        """
        Check whether the acyl carbon (acyl_c_idx) appears to be part of a peptide bond.
        Heuristically, if the acyl carbon (which bears a double bond oxygen) is also bonded to
        a nitrogen that is part of an L-alpha core (other than the candidate nitrogen), we consider it peptide-bonded.
        """
        acyl_carbon = mol.GetAtomWithIdx(acyl_c_idx)
        for nbr in acyl_carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 7 and nbr.GetIdx() != candidate_n_idx:
                # Check if this other nitrogen is near an alpha core by seeing if it is bonded to a chiral carbon
                for nn in nbr.GetNeighbors():
                    if nn.GetIdx() in alpha_indices:
                        return True
        return False

    def is_acylated_nitrogen(n_atom):
        """
        Given a nitrogen atom, return True if it is acylated.
        That is, if it has a neighboring carbon (other than the carboxyl carbon)
        which in turn bears a double-bonded oxygen.
        Also, if that acyl carbon appears part of a peptide bond (see is_peptide_bond), we discard it.
        """
        for nbr in n_atom.GetNeighbors():
            # Must be a carbon different from the carboxyl carbon.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != acid_carbon_idx:
                # Check if the carbon has a double bonded oxygen.
                for nn in nbr.GetNeighbors():
                    if nn.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                            # Exclude if the acyl carbon seems part of a peptide bond.
                            if is_peptide_bond(nbr.GetIdx(), n_atom.GetIdx()):
                                continue
                            return True
        return False

    # --- Step 4: Search for an acylated nitrogen.
    max_distance = 5  # maximum allowed bonds between the α‑carbon and the nitrogen
    found_acyl = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # nitrogen
            try:
                path = rdmolops.GetShortestPath(mol, alpha_idx, atom.GetIdx())
            except:
                continue
            if 0 < len(path)-1 <= max_distance:
                if is_acylated_nitrogen(atom):
                    found_acyl = True
                    break

    if found_acyl:
        return True, ("Molecule has an L-alpha-amino acid core (alpha carbon index %d with attached carboxyl group) "
                      "and an acylated nitrogen (within %d bonds) that is not part of a peptide bond." % (alpha_idx, max_distance))
    else:
        return False, "No qualifying N-acyl substituent found attached to the L-alpha-amino acid core."

# Example testing -- uncomment the lines below to run local tests:
# test_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"  # N-acetyl-L-aspartic acid
# result, reason = is_N_acyl_L_alpha_amino_acid(test_smiles)
# print(result, reason)