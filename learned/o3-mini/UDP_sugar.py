"""
Classifies: CHEBI:17297 UDP-sugar
"""
#!/usr/bin/env python
"""
Classifies: UDP-sugar 
Definition:
    A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to an 
    unspecified sugar via an anomeric diphosphate linkage.
    
Our strategy:
  1. Look for a UDP “core” fragment – that is, the uridine (uracil attached to a ribose)
     in which the ribose’s 5′-position bears a diphosphate group.
  2. Once the UDP core is found, require that the remainder of the molecule (i.e. 
     the part not accounted for by the UDP fragment) has at least one additional sugar ring.
     
Note: The SMARTS for the UDP core fragment is a compromise. We require that the molecule 
has a uracil ring attached to a ribose having a CH2–diphosphate substituent. The diphosphate 
is defined with two phosphorous atoms (using P with double-bonded O’s and single O’s).
"""

from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    The method requires the presence of a UDP moiety – defined as a uracil attached to a ribose
    with a 5'-diphosphate group – plus an additional sugar moiety indicated by at least one sugar ring 
    (a ring of 5 or 6 atoms containing at least one oxygen) that is not part of the UDP core.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a UDP-sugar, False otherwise
        str: A reason explaining the classification decision
    """
    
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ------------------------ 
    # 1. Search for the UDP core fragment.
    #
    # In a UDP-sugar the uridine part is attached to a diphosphate on the 5' position.
    # The following SMARTS attempts to capture the group: 
    #   – A ribose ring (5-membered) with one exocyclic CH2 group
    #   – The CH2 group carries a diphosphate comprising two phosphate units:
    #         OP(O)(=O)OP(O)(=O)O
    #   – The anomeric carbon of the ribose is attached to a uracil ring (pyrimidine with two carbonyls)
    #
    # This UDP core SMARTS does not cover all substitutions but is meant to improve our connectivity criteria.
    udp_smarts = ("O[C@H]1OC(COP(O)(=O)OP(O)(=O)O)[C@@H](O)[C@@H]1"
                  "n1ccc(=O)[nH]c1=O")
    udp_core = Chem.MolFromSmarts(udp_smarts)
    if udp_core is None:
        return False, "Error creating the UDP SMARTS fragment"
    
    udp_matches = mol.GetSubstructMatches(udp_core)
    if not udp_matches:
        return False, "UDP core fragment (uridine-diphosphate) not found"
    
    # For our present approach we take the first matching UDP fragment.
    udp_match_atoms = set(udp_matches[0])
    
    # ------------------------ 
    # 2. Search for an additional sugar ring.
    #
    # We assume that a sugar ring is typically a 5- or 6-membered ring that has at least one oxygen.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    sugar_ring_found = False
    sugar_ring_count = 0
    for ring in rings:
        # Check if ring size is 5 or 6
        if len(ring) in (5, 6):
            # Count oxygen atoms in the ring.
            oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
            if len(oxygens) >= 1:
                # We require that the sugar ring is not simply the ribose of the UDP core.
                # If not all atoms of the ring are in the UDP fragment match, then we count it as the extra sugar.
                if not set(ring).issubset(udp_match_atoms):
                    sugar_ring_found = True
                    sugar_ring_count += 1
                    
    if not sugar_ring_found:
        return False, "No additional sugar ring found outside of the UDP core fragment"
    
    # All criteria passed: UDP core is detected and an extra sugar ring exists.
    return True, ("Contains a uridine-diphosphate (UDP) core and at least one additional sugar ring, "
                  "consistent with a UDP-sugar (nucleotide attached via an anomeric diphosphate linkage)")

# For quick testing (can be removed or adapted in production)
if __name__ == "__main__":
    # Testing on a known UDP-sugar (UDP-D-glucose)
    test_smiles = ("OC[C@H]1OC(COP(O)(=O)OP(O)(=O)O)[C@H](O)[C@@H](O)[C@@H]1n1ccc(=O)[nH]c1=O")
    result, reason = is_UDP_sugar(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)