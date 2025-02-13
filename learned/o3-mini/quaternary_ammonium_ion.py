"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: Quaternary Ammonium Ion
A quaternary ammonium ion is defined as a derivative of ammonium, NH4(+),
in which all four hydrogens have been replaced by univalent (usually organyl) groups.
Thus, in our approach we look for a nitrogen (atomic number 7) that
• Has a formal charge of +1.
• Is connected to exactly four atoms with no (explicit or implicit) hydrogens.
• Has sp3 hybridization.
Then, to improve specificity we:
  - Discard candidate nitrogens that are completely embedded inside rings (i.e. all substituents are in rings).
  - In larger molecules we require that the candidate be “choline‐like” – i.e. nearby (within 4 bonds) a phosphate (P) atom.
  
If any candidate nitrogen passes, we return True along with a message including its atom index.
If no candidate is found, we return False along with the reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a quaternary ammonium ion based on its SMILES string,
    with additional filters to avoid false positives in large, rigid molecules.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule contains a (free or choline-like) quaternary ammonium ion,
              False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pre-calculate some helpful information:
    # List of phosphorus atom indices (used to check for nearby phosphate groups)
    p_atom_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
    # Total number of atoms (to help judge “small” vs “large” molecules)
    num_atoms = mol.GetNumAtoms()
    
    # Define a helper to compute topological distance (number of bonds)
    def distance_between(a_idx, b_idx):
        # Use RDKit’s shortest path function, which returns a list of atom indices in the path.
        paths = Chem.GetShortestPath(mol, a_idx, b_idx)
        if paths:
            return len(paths) - 1  # bonds count
        else:
            return None

    # Iterate over candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        # Check for formal +1 charge.
        if atom.GetFormalCharge() != 1:
            continue
        # Check that degree is 4.
        if atom.GetDegree() != 4:
            continue
        # Check that there are no attached hydrogens (explicit or implicit)
        if atom.GetTotalNumHs() != 0:
            continue
        # Check that the nitrogen is sp3 hybridized (typical for quaternary ammonium)
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Now check that not all substituents are in rings.
        neighbors = atom.GetNeighbors()
        # If every neighbor is in a ring, skip this candidate.
        if all(neighbor.IsInRing() for neighbor in neighbors):
            continue

        # For larger molecules, try to see if the candidate N is “choline‐like”
        # That is, if a phosphorus atom is present in the molecule then one expects that the quaternary N
        # (as in phosphatidylcholines) will be within a short topological path of a P atom.
        if p_atom_indices and num_atoms > 35:
            close_to_p = False
            for p_idx in p_atom_indices:
                dist = distance_between(atom.GetIdx(), p_idx)
                if dist is not None and dist <= 4:
                    close_to_p = True
                    break
            if not close_to_p:
                # Not near a phosphate group; in a large molecule this candidate is less likely to be the desired ion.
                continue
        
        # If we reached here, we found a candidate that passes our tests.
        return True, f"Found quaternary ammonium ion at atom index {atom.GetIdx()}"

    # If no candidate passed, return False.
    return False, "No suitable quaternary ammonium ion found in the molecule"

# Testing (if desired) using one of the provided examples:
if __name__ == "__main__":
    # Example: bretylium (expected: True)
    test_smiles = "CC[N+](C)(C)Cc1ccccc1Br"
    result, reason = is_quaternary_ammonium_ion(test_smiles)
    print(result, reason)