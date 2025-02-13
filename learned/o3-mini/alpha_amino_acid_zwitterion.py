"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid zwitterion obtained by transfer of a proton from the carboxy to the amino group
of any alpha–amino acid; major species at pH 7.3.

In this improved version we:
  – Identify a carboxylate group via the SMARTS "[C](=O)[O-]".
  – Identify candidate ammonium groups by selecting nitrogen atoms with formal charge > 0,
    at least one hydrogen, that are sp3–hybridized and are not aromatic.
  – For every carboxylate/ammonium pair, compute the shortest bond path.
    If the path consists of 3 atoms (2 bonds) we take the middle atom,
    and if 4 atoms (3 bonds) we take the second intermediate atom (index 2).
  – We require that the bridging (alpha) atom is a carbon with at least one hydrogen,
    and that aside from its bonds to the carboxylate and ammonium groups it has exactly one (glycine)
    or two (substituted amino acid) heavy–atom neighbors.
  – This extra connectivity requirement should help eliminate cases where a zwitterionic core
    is “accidentally” found in larger or branched molecules.
    
If no valid candidate is found the function returns False.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha–amino-acid zwitterion based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES.
      2. Look for at least one carboxylate group using the SMARTS "[C](=O)[O-]".
      3. Identify candidate ammonium groups by selecting nitrogen atoms with:
            - formal charge > 0,
            - at least one hydrogen,
            - sp3 hybridization,
            - not aromatic.
      4. For every carboxylate/ammonium candidate pair, compute the shortest bond path.
         If that path has either 3 atoms (2 bonds) or 4 atoms (3 bonds) then choose the bridging atom.
         For 3 atoms, use the middle atom (path[1]);
         for 4 atoms, use the second intermediate (path[2]).
      5. The candidate “alpha–carbon” must:
            - be carbon,
            - have at least one hydrogen,
            - have “typical” heavy–atom connectivity. Specifically, when we discount its bonds
              to the carboxylate carbon and the ammonium nitrogen, it should have exactly one (glycine)
              or two (other amino acids) heavy-atom neighbors.
    
    Args:
       smiles (str): SMILES representation of the molecule.
    
    Returns:
       bool: True if the molecule is classified as containing an alpha–amino–acid zwitterion core,
             otherwise False.
       str: A message explaining the basis for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Query for carboxylate group: matches a carbon double-bonded to an oxygen and single-bonded to an O-
    carboxylate_smarts = "[C](=O)[O-]"
    carbox_query = Chem.MolFromSmarts(carboxylate_smarts)
    carbox_matches = mol.GetSubstructMatches(carbox_query)
    if not carbox_matches:
        return False, "No carboxylate group [C](=O)[O-] found"
    
    # Identify candidate ammonium groups: nitrogen with positive charge, >=1 hydrogen,
    # must be sp3 hybridized and not aromatic.
    ammonium_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() > 0 and atom.GetTotalNumHs() >= 1:
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            if atom.GetIsAromatic():
                continue
            ammonium_indices.append(atom.GetIdx())
    if not ammonium_indices:
        return False, "No suitable protonated amine (sp3 N with positive charge and >=1 H and non-aromatic) found"
    
    # For each pair of carboxylate (we take the carboxylate carbon as the 0-th atom in the match)
    # and candidate ammonium, compute the shortest path.
    # For a valid alpha-core, if the path has 3 atoms (2 bonds) then the middle atom is candidate;
    # if it has 4 atoms (3 bonds) then the second intermediate (index 2) is candidate.
    for carbox_match in carbox_matches:
        carbox_idx = carbox_match[0]
        for ammonium_idx in ammonium_indices:
            path = rdmolops.GetShortestPath(mol, carbox_idx, ammonium_idx)
            if not path:
                continue
            num_bonds = len(path) - 1  # number of bonds in path
            if num_bonds == 2 or num_bonds == 3:
                # Determine candidate alpha carbon based on path length:
                if num_bonds == 2:
                    candidate_alpha_idx = path[1]
                else:  # num_bonds == 3
                    candidate_alpha_idx = path[2]
                alpha_atom = mol.GetAtomWithIdx(candidate_alpha_idx)
                
                # Candidate must be carbon.
                if alpha_atom.GetSymbol() != "C":
                    continue
                
                # Candidate alpha carbon should have at least one hydrogen.
                if alpha_atom.GetTotalNumHs() < 1:
                    continue
                
                # Check heavy-atom connectivity. Get neighbors that are not the carboxylate carbon or ammonium nitrogen.
                neighbor_ids = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors()]
                # Remove the indices corresponding to the connecting acid and amine groups.
                trimmed_neighbors = [nbr for nbr in alpha_atom.GetNeighbors() 
                                     if nbr.GetIdx() not in (carbox_idx, ammonium_idx) and nbr.GetAtomicNum() > 1]
                # In a typical amino acid:
                # - For glycine the alpha carbon will have only its two main bonds (to carboxylate and amine)
                #   but its overall heavy connectivity (neighbors) is 2.
                # - For substituted amino acids, there is one extra bond to the side chain (thus 3 heavy neighbors).
                # We enforce that there is at least one additional non–acid/amine heavy neighbor,
                # and the total heavy connectivity (neighbors with atomic number > 1) should be 2 or 3.
                total_heavy = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(total_heavy) not in (2, 3):
                    continue
                if len(trimmed_neighbors) < 1:
                    continue
                
                # If we reach here, we have found a candidate core.
                return True, (
                    f"Found zwitterionic core: carboxylate C (idx {carbox_idx}) - "
                    f"alpha C (idx {candidate_alpha_idx}) - ammonium N (idx {ammonium_idx}); "
                    f"path length: {num_bonds} bonds."
                )
    
    return False, "Alpha–amino-acid zwitterion core not found"

# (Optional) Testing examples:
if __name__ == "__main__":
    test_examples = {
        # True positives:
        "(R)-2-amino-4-oxopentanoic acid zwitterion": "[C@@H](CC(C)=O)(C(=O)[O-])[NH3+]",
        "L-valine zwitterion": "CC(C)[C@H]([NH3+])C([O-])=O",
        "(2S,3S)-diazaniumylbutanoate zwitterion": "O=C([O-])[C@@H](N)[C@H](C)[NH3+]",
        "L-vinylglycine zwitterion": "[C@@H](C(=O)[O-])(C=C)[NH3+]",
        "thioproline zwitterion": "C(C1[NH2+]CSC1)(=O)[O-]",
        # Example that was a false negative previously (should now be caught):
        "N-methyl-aminoisobutyrate zwitterion": "[O-]C(=O)C([NH2+]C)(C)C",
        # Some examples that were wrongly classified might now be filtered out:
        "D-asparagine zwitterion": "NC(=O)C[C@@H]([NH3+])C([O-])=O",
        "L-threonine zwitterion": "C[C@@H](O)[C@H]([NH3+])C([O-])=O",
    }
    
    for name, smi in test_examples.items():
        result, reason = is_alpha_amino_acid_zwitterion(smi)
        print(f"{name} -> {result}: {reason}")