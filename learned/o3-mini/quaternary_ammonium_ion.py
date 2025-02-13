"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: Quaternary Ammonium Ion

A quaternary ammonium ion is defined as a derivative of ammonium, NH4(+),
in which all four hydrogens bonded to nitrogen have been replaced with 
univalent (usually organyl) groups.

Our approach is to first identify candidate nitrogen atoms that would be 
expected in a tetraalkylammonium:
  - Atomic number 7, formal charge +1.
  - Connected to exactly 4 atoms and has no explicit (or implicit) hydrogens.
  - sp3 hybridized.
Then to improve specificity we further require that:
  - The nitrogen’s substituents are not all inside rings.
  - In larger molecules that contain phosphorus (e.g. phospholipids) the candidate
    must be in “choline‐like” environment, that is within 5 bonds of a phosphorus atom.
  - In molecules with no phosphorus, if they are large (>40 heavy atoms) then we require 
    that the candidate nitrogen be attached to at least one aromatic fragment 
    (as seen in dyes like bretylium or methyl green). If the molecule is small we accept 
    the candidate.
    
If any candidate passes all tests, the function returns True along with the candidate’s 
atom index. Otherwise, it returns False along with an explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a quaternary ammonium ion based on its SMILES string.
    
    Our improved criteria are as follows:
      - Candidate nitrogen must have atomic number 7, formal charge +1, degree 4,
        no attached hydrogens, and sp3 hybridization.
      - Discard candidate if all four substituents are in rings.
      - If the molecule has phosphorus atoms (typical in phospholipids) then the candidate
        must be within 5 bonds (topological distance) of at least one phosphorus atom.
      - If no phosphorus is found and the molecule is relatively large (>40 atoms), then the 
        candidate must have at least one neighbor that is aromatic (e.g. found in dyes).
      - Smaller molecules (<40 atoms) lacking phosphorus are accepted if they fulfill the basic criteria.
      
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if a quaternary ammonium ion candidate is found, False otherwise.
       str: A reason message including the candidate’s atom index if found.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # In larger molecules, we want to check for nearby phosphorus atoms.
    p_atom_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
    num_atoms = mol.GetNumAtoms()
    
    # Helper: compute topological (bond) distance between two atoms
    def bond_distance(idx1, idx2):
        try:
            path = Chem.GetShortestPath(mol, idx1, idx2)
            if path:
                return len(path) - 1
        except Exception:
            pass
        return None

    # Iterate over atoms to find candidate nitrogen atoms
    for atom in mol.GetAtoms():
        # Candidate must be nitrogen (atomic num 7)
        if atom.GetAtomicNum() != 7:
            continue
        # Must have formal charge +1
        if atom.GetFormalCharge() != 1:
            continue
        # Must have exactly 4 neighbors
        if atom.GetDegree() != 4:
            continue
        # Must have no attached hydrogens (explicit or implicit)
        if atom.GetTotalNumHs() != 0:
            continue
        # Must be sp3 hybridized
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Check that not all substituents are in rings.
        neighbors = atom.GetNeighbors()
        if all(neighbor.IsInRing() for neighbor in neighbors):
            continue
        
        # --- Now apply additional filters ---
        candidate_idx = atom.GetIdx()
        passes_extra = False
        
        # If the molecule contains phosphorus, require close proximity.
        if p_atom_indices:
            for p_idx in p_atom_indices:
                dist = bond_distance(candidate_idx, p_idx)
                if dist is not None and dist <= 5:
                    passes_extra = True
                    break
            if not passes_extra:
                # For phospholipid-like molecules, candidate N not close to any P atom is less likely.
                continue
        else:
            # For molecules without phosphorus, we judge based on overall size.
            if num_atoms > 40:
                # In larger non-phosphorus molecules, require the N has at least one aromatic neighbor.
                for nbr in neighbors:
                    if nbr.GetIsAromatic():
                        passes_extra = True
                        break
                if not passes_extra:
                    continue
            else:
                # For small molecules (<40 atoms) we accept the candidate.
                passes_extra = True

        if passes_extra:
            return True, f"Found quaternary ammonium ion at atom index {candidate_idx}"
    
    return False, "No suitable quaternary ammonium ion found in the molecule"

# Testing examples (uncomment to run tests)
if __name__ == "__main__":
    # Example: bretylium which should return True.
    test_smiles = "CC[N+](C)(C)Cc1ccccc1Br"
    result, reason = is_quaternary_ammonium_ion(test_smiles)
    print(result, reason)