"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: UDP-sugar
Definition:
    A pyrimidine nucleotide-sugar having UDP as the nucleotide component 
    (a uracil attached to a ribose bearing a 5'-diphosphate group) linked to an 
    unspecified sugar via an anomeric diphosphate linkage.
    
Our strategy revised:
  1. Look for a diphosphate fragment defined by the SMARTS "OP(O)(=O)OP(O)(=O)O".
  2. Look for a uracil fragment defined by the SMARTS "n1ccc(=O)[nH]c1=O".
  3. Verify that a diphosphate atom and a uracil atom are connected by a short bond path.
  4. Then, require that in addition to the UDP core (which includes the ribose),
     there is at least one extra sugar ring (a 5- or 6-membered ring that contains oxygen)
     outside of the UDP core.
     
We use the revised uracil SMARTS to cover alternate representations.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    
    The method requires the presence of a UDP moiety – defined by:
        a) a diphosphate fragment (OP(O)(=O)OP(O)(=O)O)
        b) a uracil fragment (n1ccc(=O)[nH]c1=O)
    that are connected through a short bond‐path (≤6 bonds, taken as a proxy for the ribose)
    and that the molecule contains at least one additional sugar ring (5- or 6-membered ring with at least one oxygen)
    outside of the UDP core.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a UDP-sugar; False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---------- Step 1: Identify diphosphate fragment ----------
    # Use SMARTS for a diphosphate fragment.
    dp_smarts = "OP(O)(=O)OP(O)(=O)O"
    dp_mol = Chem.MolFromSmarts(dp_smarts)
    if dp_mol is None:
        return False, "Error creating diphosphate SMARTS fragment"
    dp_matches = mol.GetSubstructMatches(dp_mol, useChirality=False)
    if not dp_matches:
        return False, "Diphosphate fragment not found"
    
    # ---------- Step 2: Identify uracil fragment ----------
    # Revised SMARTS pattern for the uracil moiety.
    uracil_smarts = "n1ccc(=O)[nH]c1=O"
    uracil_mol = Chem.MolFromSmarts(uracil_smarts)
    if uracil_mol is None:
        return False, "Error creating uracil SMARTS fragment"
    uracil_matches = mol.GetSubstructMatches(uracil_mol, useChirality=False)
    if not uracil_matches:
        return False, "Uracil fragment not found"
    
    # ---------- Step 3: Verify UDP core connectivity (diphosphate linked to uracil via a short path) ----------
    found_udp_core = False
    udp_core_atoms = set()
    cutoff = 6  # maximum bond distance (in bonds)
    for dp_match in dp_matches:
        for ur_match in uracil_matches:
            # Check each pair of atoms from dp_match and ur_match.
            for dp_atom in dp_match:
                for ur_atom in ur_match:
                    try:
                        path = rdmolops.GetShortestPath(mol, dp_atom, ur_atom)
                    except Exception:
                        continue
                    if path and (len(path) - 1) <= cutoff:
                        udp_core_atoms.update(dp_match)
                        udp_core_atoms.update(ur_match)
                        udp_core_atoms.update(path)
                        found_udp_core = True
                        break
                if found_udp_core:
                    break
            if found_udp_core:
                break
        if found_udp_core:
            break
    if not found_udp_core:
        return False, "UDP core fragment (uridine-diphosphate) not found or not properly connected"
    
    # ---------- Step 4: Look for an additional sugar ring outside of the UDP core ----------
    # A sugar ring is defined as a 5- or 6-membered ring containing at least one oxygen.
    ring_info = mol.GetRingInfo()
    sugar_ring_found = False
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            # Count oxygen atoms in the ring.
            oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxygen_count:
                # Ensure ring is not completely contained within UDP core.
                if not set(ring).issubset(udp_core_atoms):
                    sugar_ring_found = True
                    break
    if not sugar_ring_found:
        return False, "No additional sugar ring found outside of the UDP core"
    
    return True, ("Contains a UDP core (diphosphate linked to a uracil-bearing ribose) and an extra sugar ring, "
                  "consistent with a UDP-sugar structure.")

# Quick test (this can be removed or modified as needed)
if __name__ == "__main__":
    test_smiles = "OC[C@H]1OC(COP(O)(=O)OP(O)(=O)O)[C@H](O)[C@@H](O)[C@@H]1n1ccc(=O)[nH]c1=O"
    result, reason = is_UDP_sugar(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)