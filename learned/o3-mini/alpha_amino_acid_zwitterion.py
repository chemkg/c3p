"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid-zwitterion obtained by transfer of a proton from the carboxy to the amino group 
of any alpha-amino acid; major species at pH 7.3.
In an idealized alpha amino acid zwitterion the carboxylate (–C(=O)[O–]) and the ammonium ([NH3+])
originate from a neutral amino acid (NH2–CHR–COOH) so that the core is R–C(H)(NH3+)(C(=O)[O–]).
This program attempts to detect such a motif by:
  (1) Finding carboxylate groups (via SMARTS "[C](=O)[O-]") and ammonium groups (via SMARTS "[NH3+]").
  (2) For each carboxylate/ammonium pair, computing the shortest path.
  (3) Requiring that the carboxylate carbon and ammonium nitrogen are connected by either 2 bonds
      (the canonical case, where both groups are directly attached to the same alpha carbon) or by 3 bonds
      (e.g. in some di‐zwitterions) and that the intervening alpha carbon has at least one hydrogen.
This additional connectivity check is intended to reduce false positives (molecules that have both an –COO– and [NH3+]
but not arranged as a true amino acid) and false negatives.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    It first checks that the molecule contains a carboxylate group (C(=O)[O-])
    and a protonated ammonium (NH3+). Then it searches for a candidate “alpha‐carbon”
    that bridges the carboxylate and ammonium groups. For a classical alpha-amino acid zwitterion,
    the carboxylate carbon and the ammonium nitrogen will be connected by exactly two bonds
    (i.e. via a central alpha-carbon) or in some cases three bonds.
    
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
    
    # First we check that the molecule contains a carboxylate group.
    # This SMARTS identifies a carbon that is double-bonded to an oxygen and single-bonded to an O–.
    carboxylate_smarts = "[C](=O)[O-]"
    carboxylate_query = Chem.MolFromSmarts(carboxylate_smarts)
    carbox_matches = mol.GetSubstructMatches(carboxylate_query)
    if not carbox_matches:
        return False, "No carboxylate group [C](=O)[O-] found"
    
    # Next check for a protonated amino group.
    ammonium_smarts = "[NH3+]"
    ammonium_query = Chem.MolFromSmarts(ammonium_smarts)
    ammonium_matches = mol.GetSubstructMatches(ammonium_query)
    if not ammonium_matches:
        return False, "No protonated ammonium group [NH3+] found"
    
    # For each found carboxylate group, use its carbon atom as candidate donor.
    # For each found ammonium, compute the shortest bond path between the carboxylate carbon and ammonium nitrogen.
    for carbox_match in carbox_matches:
        carbox_idx = carbox_match[0]  # the carbon index in the carboxylate
        for ammonium_match in ammonium_matches:
            ammonium_idx = ammonium_match[0]
            # Compute the shortest path (list of atom indices) between the two atoms.
            path = rdmolops.GetShortestPath(mol, carbox_idx, ammonium_idx)
            # For a proper amino acid, the carboxylate and ammonium should be connected through the alpha-carbon.
            # In the canonical case the path should be: carboxylate C -- alpha C -- ammonium N (3 atoms, 2 bonds).
            # In some zwitterions a one-atom spacer might be inserted so that the number of bonds is 3.
            num_bonds = len(path) - 1
            if num_bonds == 2 or num_bonds == 3:
                # Identify the candidate alpha carbon.
                # If num_bonds==2 then the atom common to both groups is the alpha carbon;
                # if num_bonds==3 then the middle atom (path[1]) is our candidate.
                if num_bonds == 2:
                    alpha_idx = path[1]
                else:
                    alpha_idx = path[1]  # for a 3-bond path, we take the first intermediate atom as the candidate
                
                alpha_atom = mol.GetAtomWithIdx(alpha_idx)
                # Check that the candidate alpha atom is indeed a carbon
                if alpha_atom.GetSymbol() != "C":
                    continue
                # And that it has at least one hydrogen (often alpha carbons contain one H).
                if alpha_atom.GetTotalNumHs() < 1:
                    continue
                    
                # If we reach here, we have a candidate alpha carbon bridging a carboxylate and ammonium.
                # Return True with a message including the bond path found.
                return True, f"Found zwitterionic core: carboxylate C (idx {carbox_idx}) - alpha C (idx {alpha_idx}) - ammonium N (idx {ammonium_idx}); path length: {num_bonds} bonds."
    
    # If no bridged connection is found, the molecule fails the classifier.
    return False, "Alpha-amino-acid zwitterion core not found"


# (Optional) Testing examples:
if __name__ == "__main__":
    # Examples from the provided list:
    examples = {
        "(R)-2-amino-4-oxopentanoic acid zwitterion": "[C@@H](CC(C)=O)(C(=O)[O-])[NH3+]",
        "L-valine zwitterion": "CC(C)[C@H]([NH3+])C([O-])=O",
        "(2S,3S)-diazaniumylbutanoate zwitterion": "O=C([O-])[C@@H](N)[C@H](C)[NH3+]",
        "N-methyl-aminoisobutyrate zwitterion": "[O-]C(=O)C([NH2+]C)(C)C"
    }
    
    for name, smi in examples.items():
        result, reason = is_alpha_amino_acid_zwitterion(smi)
        print(f"{name}: {result} ({reason})")