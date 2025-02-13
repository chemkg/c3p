"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
A valid epoxy fatty acid must:
  - Have a carboxylic acid group (fatty acid functionality)
  - Contain a long, largely linear aliphatic chain (here we check that the chain 
    attached at the acid end has a sufficient number of allowed carbon atoms)
  - Contain at least one epoxide ring, defined as a three‐membered cycle [C;r3][O;r3][C;r3]
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    
    The strategy is as follows:
      1. Parse the SMILES.
      2. Check that a carboxylic acid group is present.
      3. Check that the overall molecule has a sufficient number of carbon atoms.
      4. Look for a three‐membered epoxide ring.
      5. Starting at the carboxylic acid carbon, do a depth–first search 
         (only along carbons that are “allowed” – i.e. not in rings larger than 3 atoms)
         to compute the longest linear chain. This approximates the fatty acid chain.
      6. If these criteria are met then return True.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as an epoxy fatty acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (1) Check for the carboxylic acid group.
    # Use a SMARTS that catches a typical free -COOH unit.
    acid_smarts = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    if not acid_matches:
        return False, "Missing carboxylic acid group (fatty acid functionality)"
    
    # (2) Check overall number of carbon atoms.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 18:
        return False, f"Too few carbon atoms ({len(c_atoms)}); not long enough to be a fatty acid"
    
    # (3) Identify at least one epoxide ring.
    epoxide_pattern = Chem.MolFromSmarts("[C;r3][O;r3][C;r3]")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring (three-membered heterocycle containing oxygen) detected."
        
    # (4) Check for a long linear chain emanating from the acid group.
    # For fatty acids the carboxylic acid is typically at one terminus of a long aliphatic chain.
    # We define a helper: allowed chain atoms are carbons that are either not in any ring
    # or if they are in a ring, the ring size must be exactly 3 (allowing for the epoxide).
    ring_info = mol.GetRingInfo()
    # Precompute a list of rings (as sets of atom indices) for efficiency.
    rings = [set(r) for r in ring_info.AtomRings()]
    
    def in_large_ring(atom_idx):
        # Returns True if the atom is in any ring larger than 3 atoms.
        for r in rings:
            if atom_idx in r and len(r) > 3:
                return True
        return False

    def is_allowed_chain_atom(atom):
        # We only allow carbon atoms for the chain.
        # If the atom is in a ring, then it is allowed only if it is exclusively
        # in three-membered rings (if any). (A pragmatic approximation.)
        if atom.GetAtomicNum() != 6:
            return False
        # If not in any ring, allow.
        atom_idx = atom.GetIdx()
        if not mol.GetAtomWithIdx(atom_idx).IsInRing():
            return True
        # If in a ring, check all rings that contain the atom:
        for r in rings:
            if atom_idx in r and len(r) > 3:
                return False
        return True

    # We now perform a DFS starting from the carbon atom that is the "acid carbon".
    # (The SMARTS "C(=O)[O;H]" returns a tuple with [acid carbon, oxygen].)
    # We will follow bonds only through atoms that meet is_allowed_chain_atom.
    def dfs(atom_idx, visited):
        max_length = 0
        current_atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in current_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            if is_allowed_chain_atom(nbr):
                visited.add(nbr_idx)
                path_length = 1 + dfs(nbr_idx, visited)
                if path_length > max_length:
                    max_length = path_length
                visited.remove(nbr_idx)
        return max_length

    longest_chain = 0
    # Try each acid group match and consider the acid carbon (first index in the match).
    for match in acid_matches:
        acid_carbon_idx = match[0]
        acid_atom = mol.GetAtomWithIdx(acid_carbon_idx)
        # Look at neighbors of the acid carbon that are allowed chain atoms.
        for nbr in acid_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and is_allowed_chain_atom(nbr):
                visited = {acid_carbon_idx, nbr.GetIdx()}
                chain_length = 1 + dfs(nbr.GetIdx(), visited)
                if chain_length > longest_chain:
                    longest_chain = chain_length

    # Require that the chain (not including the acid carbon itself) is long enough.
    # This threshold can be tuned – here we require at least 12 consecutive carbons.
    if longest_chain < 12:
        return False, f"Aliphatic chain too short (found chain length {longest_chain}, need at least 12 carbons)"
    
    return True, "Molecule has a carboxylic acid group, a sufficiently long aliphatic chain, and an epoxide ring."

# Testing examples if run as main:
if __name__ == "__main__":
    # One of the true positives:
    test_smiles = "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O"
    result, reason = is_epoxy_fatty_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)