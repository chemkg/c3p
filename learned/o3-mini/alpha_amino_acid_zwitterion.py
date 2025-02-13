"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid zwitterion obtained by transfer of a proton from the carboxyl to the amino group
of any alpha-amino acid; major species at pH 7.3.

This program uses two main substructure queries:
  - a carboxylate group – identified by the SMARTS "[C](=O)[O-]"
  - an ammonium group – now defined more broadly as any nitrogen atom with a positive formal charge
    and at least one hydrogen (to capture [NH3+], [NH2+], etc.).
For every carboxylate/ammonium pair found in the molecule we calculate the shortest bond path.
If that path has either 2 bonds (3 atoms in the path) or 3 bonds (4 atoms in the path) then the
atom bridging the two charged groups is taken as the candidate “alpha–carbon.”
For an idealized alpha–amino acid zwitterion this intermediate atom should be a carbon that (a) carries
at least one hydrogen and (b) has the expected heavy-atom connectivity (typically 2 neighbors for glycine
or 3 for substituted amino acids). These extra criteria were added in an attempt to eliminate cases
where the two charged groups occur elsewhere and thereby reduce false classifications.
If such a candidate core is identified the function returns True along with details of the match.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES.
      2. Look for at least one carboxylate group using the SMARTS "[C](=O)[O-]".
      3. Identify any nitrogen atoms that are protonated – i.e. those with formal charge > 0 and with >= 1 hydrogen.
      4. For every carboxylate and ammonium candidate pair, compute the shortest bond path.
         We accept a candidate if the number of bonds in the path is 2 or 3.
      5. In that path the intermediate atom is the candidate alpha–carbon. It must be a carbon and should have
         at least one attached hydrogen plus “typical” connectivity (i.e. 2 or 3 heavy neighbors).
    
    Args:
       smiles (str): SMILES representation of the molecule.
    
    Returns:
       bool: True if the molecule is classified as containing an alpha-amino-acid zwitterion core,
             otherwise False.
       str: A message explaining the basis for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Query for carboxylate group: carbon bonded to a carbonyl oxygen and an anionic oxygen.
    carboxylate_smarts = "[C](=O)[O-]"
    carbox_query = Chem.MolFromSmarts(carboxylate_smarts)
    carbox_matches = mol.GetSubstructMatches(carbox_query)
    if not carbox_matches:
        return False, "No carboxylate group [C](=O)[O-] found"
    
    # Identify candidate ammonium groups.
    ammonium_indices = []
    for atom in mol.GetAtoms():
        # Select nitrogen atoms with positive charge and with at least one hydrogen
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() > 0 and atom.GetTotalNumHs() >= 1:
            ammonium_indices.append(atom.GetIdx())
    if not ammonium_indices:
        return False, "No protonated amine (N with positive charge and at least one hydrogen) found"
    
    # Now loop over all pairs (carboxylate, ammonium) and check the connectivity.
    # Each carboxylate match: we assume the first atom (index 0 in the tuple) is the carboxylate carbon.
    for carbox_match in carbox_matches:
        carbox_idx = carbox_match[0]
        for ammonium_idx in ammonium_indices:
            # Compute the shortest bond path between the carboxylate carbon and the ammonium nitrogen.
            path = rdmolops.GetShortestPath(mol, carbox_idx, ammonium_idx)
            if not path:
                continue  # no connection found
            num_bonds = len(path) - 1
            # Accept cases where path has 2 or 3 bonds only.
            if num_bonds == 2 or num_bonds == 3:
                # For both cases, take the first intermediate atom as the candidate alpha–carbon.
                # For a 2-bond (3-atom) path, path[1] is directly between the two.
                candidate_alpha_idx = path[1]
                alpha_atom = mol.GetAtomWithIdx(candidate_alpha_idx)
                # The candidate must be a carbon.
                if alpha_atom.GetSymbol() != "C":
                    continue
                # Must have at least one hydrogen.
                if alpha_atom.GetTotalNumHs() < 1:
                    continue
                # Check the heavy atom connectivity: in a typical alpha amino acid,
                # the alpha–carbon is attached to either 2 (glycine case) or 3 heavy atoms.
                heavy_neighbors = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(heavy_neighbors) not in (2, 3):
                    continue
                # We have a candidate core. Provide details.
                return True, (
                    f"Found zwitterionic core: carboxylate C (idx {carbox_idx}) - "
                    f"alpha C (idx {candidate_alpha_idx}) - ammonium N (idx {ammonium_idx}); "
                    f"path length: {num_bonds} bonds."
                )
    
    # If none of the pairs provided a valid alpha–carbon bridging the groups, then this is not a match.
    return False, "Alpha–amino-acid zwitterion core not found"

# (Optional) Testing examples:
if __name__ == "__main__":
    # A few examples from the provided list. (More comprehensive testing is needed to balance false positives/negatives.)
    examples = {
        "(R)-2-amino-4-oxopentanoic acid zwitterion": "[C@@H](CC(C)=O)(C(=O)[O-])[NH3+]",
        "L-valine zwitterion": "CC(C)[C@H]([NH3+])C([O-])=O",
        "(2S,3S)-diazaniumylbutanoate zwitterion": "O=C([O-])[C@@H](N)[C@H](C)[NH3+]",
        "N-methyl-aminoisobutyrate zwitterion": "[O-]C(=O)C([NH2+]C)(C)C",
        "thioproline zwitterion": "C(C1[NH2+]CSC1)(=O)[O-]",
        "N-methylproline zwitterion": "[C@H]1(CCC[NH+]1C)C([O-])=O"
    }
    
    for name, smi in examples.items():
        result, reason = is_alpha_amino_acid_zwitterion(smi)
        print(f"{name}: {result} ({reason})")