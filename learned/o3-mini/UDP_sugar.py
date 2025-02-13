"""
Classifies: CHEBI:17297 UDP-sugar
"""
#!/usr/bin/env python
"""
Classifies: UDP-sugar
Definition:
    A pyrimidine nucleotide-sugar having UDP as the nucleotide component 
    (a uracil attached to a ribose bearing a 5'-diphosphate group) linked to an 
    unspecified sugar via an anomeric diphosphate linkage.
    
Our strategy:
  1. Look for a diphosphate fragment defined by the SMARTS "OP(O)(=O)OP(O)(=O)O".
  2. Look for a uracil fragment defined by the SMARTS "n1c(=O)[nH]c(=O)n1".
  3. Verify that a diphosphate atom and a uracil atom are connected by a short bond path.
  4. Then, require that in addition to the UDP core (which includes the ribose),
     there is at least one extra sugar ring (a 5- or 6-membered ring that contains oxygen)
     outside of the UDP core.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    
    The method requires the presence of a UDP moiety – defined here by:
        a) a diphosphate fragment (OP(O)(=O)OP(O)(=O)O)
        b) a uracil fragment (n1c(=O)[nH]c(=O)n1)
    that are connected through a short bond‐path (≤6 bonds, taken as a proxy for the ribose)
    and that the molecule contains at least one additional sugar ring (5- or 6-membered ring with an oxygen)
    outside of the UDP core.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a UDP-sugar; False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---------- Step 1: Identify diphosphate fragment ----------
    # SMARTS for a diphosphate fragment (note: many UDP sugars show exactly this connectivity)
    dp_smarts = "OP(O)(=O)OP(O)(=O)O"
    dp_mol = Chem.MolFromSmarts(dp_smarts)
    if dp_mol is None:
        return False, "Error creating diphosphate SMARTS fragment"
    dp_matches = mol.GetSubstructMatches(dp_mol, useChirality=False)
    if not dp_matches:
        return False, "Diphosphate fragment not found"
    
    # ---------- Step 2: Identify uracil fragment ----------
    # SMARTS for the uracil moiety, without enforcing chirality.
    uracil_smarts = "n1c(=O)[nH]c(=O)n1"
    uracil_mol = Chem.MolFromSmarts(uracil_smarts)
    if uracil_mol is None:
        return False, "Error creating uracil SMARTS fragment"
    uracil_matches = mol.GetSubstructMatches(uracil_mol, useChirality=False)
    if not uracil_matches:
        return False, "Uracil fragment not found"
    
    # ---------- Step 3: Verify connectivity between diphosphate and uracil as UDP core ----------
    # We want at least one diphosphate atom and one uracil atom to be connected by a short path.
    # (This is our proxy for the ribose linking them.)
    found_udp_core = False
    udp_core_atoms = set()
    cutoff = 6  # maximum bond distance allowed between a dp atom and a uracil atom
    for dp_match in dp_matches:
        for ur_match in uracil_matches:
            # Check each pair of atoms from dp and uracil matches
            for dp_atom in dp_match:
                for ur_atom in ur_match:
                    path = rdmolops.GetShortestPath(mol, dp_atom, ur_atom)
                    if path and (len(path) - 1) <= cutoff:
                        # Build a set of atoms that comprise the UDP core:
                        # (the dp match, uracil match, and the atoms along the connecting path)
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
        return False, "UDP core fragment (uridine-diphosphate) not found or not connected"
    
    # ---------- Step 4: Look for an additional sugar ring not part of the UDP core ----------
    # Define a sugar ring as a 5- or 6-membered ring that contains at least one oxygen.
    ring_info = mol.GetRingInfo()
    sugar_ring_found = False
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            # Count oxygen atoms in the ring.
            oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
            if oxygens:
                # Ensure that this ring is not completely contained in the UDP core.
                if not set(ring).issubset(udp_core_atoms):
                    sugar_ring_found = True
                    break
                    
    if not sugar_ring_found:
        return False, "No additional sugar ring found outside of the UDP core"
    
    return True, ("Contains a UDP core (diphosphate linked to a uracil-bearing ribose) and an extra sugar ring, "
                  "consistent with a UDP-sugar structure.")

# For quick testing on one example (this section can be removed or modified for production)
if __name__ == "__main__":
    test_smiles = "OC[C@H]1OC(COP(O)(=O)OP(O)(=O)O)[C@H](O)[C@@H](O)[C@@H]1n1c(=O)[nH]c(=O)n1"
    result, reason = is_UDP_sugar(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)