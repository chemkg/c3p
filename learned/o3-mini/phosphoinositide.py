"""
Classifies: CHEBI:18179 phosphoinositide
"""
#!/usr/bin/env python
"""
Classifies: Phosphoinositide

Definition:
  A phosphoinositide is any phosphatidylinositol that is phosphorylated at one or more of 
  the hydroxy groups of the inositol ring. In addition to the extra phosphorylation on the ring,
  a phosphoinositide must present a bridging phosphate that connects the inositol ring to a diacylglycerol fragment.
  
Revised strategy:
  1. Parse the SMILES and add explicit hydrogens.
  2. Identify candidate inositol rings: six-membered rings composed entirely 
     of sp³ carbon atoms. Each ring carbon is checked for “OH-like” substituents 
     (an oxygen with either a hydrogen or bound to P). At least 4 of 6 must have such a substituent.
  3. For the candidate ring, collect phosphate (P) atoms attached through oxygen substituents.
  4. For each phosphate attached to the ring, check for a “bridging” oxygen:
     One of its non-ring oxygen neighbors must be attached to a carbonyl group [CX3](=O).
  5. Require that at least one phosphate shows the bridging feature and that 
     at least 2 distinct phosphate groups are attached.
     
Note:
  This algorithm is heuristic and may not capture every edge case.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    
    A phosphoinositide is defined as a phosphatidylinositol carrying at least one extra phosphate
    on the inositol ring, with at least one phosphate linking the inositol to a diacylglycerol fragment.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phosphoinositide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so OH groups are visible.
    mol = Chem.AddHs(mol)
    
    # --- Step 1: Identify candidate inositol ring ---
    # We require a six-membered, all-carbon, sp3 ring with at least 4 carbons bearing an “OH-like” substituent.
    ring_info = mol.GetRingInfo()
    candidate_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        valid_ring = True
        substitution_count = 0  # how many ring carbons have an 'OH-like' substituent
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                valid_ring = False
                break
            # Look in neighbors outside the ring.
            found_substituent = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check for an OH-like group: oxygen with at least one hydrogen and formal charge zero
                    if nbr.GetFormalCharge() == 0 and any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                        found_substituent = True
                        break
                    # Alternatively, if the oxygen is connected to a phosphorus.
                    if any(n.GetAtomicNum() == 15 for n in nbr.GetNeighbors()):
                        found_substituent = True
                        break
            if found_substituent:
                substitution_count += 1
        if valid_ring and substitution_count >= 4:
            candidate_ring = ring
            break
    
    if candidate_ring is None:
        return False, "No inositol-like six-membered (cyclohexane) ring with sufficient OH-like substituents was found"
    
    # --- Step 2: Gather phosphate groups attached to the candidate ring ---
    # For each ring atom, check neighbor oxygens. If an oxygen is bonded to a phosphorus, record that phosphorus.
    phosphate_set = set()
    # Also keep a mapping: phosphate index -> associated ring oxygen indices (for bridging check)
    phosphate_ring_oxygens = {}
    
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in candidate_ring:
                continue
            if nbr.GetAtomicNum() != 8:
                continue
            # Check if this oxygen is connected to any phosphorus:
            for subnbr in nbr.GetNeighbors():
                if subnbr.GetAtomicNum() == 15:
                    phosphate_set.add(subnbr.GetIdx())
                    phosphate_ring_oxygens.setdefault(subnbr.GetIdx(), []).append(nbr.GetIdx())
    
    if len(phosphate_set) < 2:
        return False, f"Only {len(phosphate_set)} phosphate group(s) found on the inositol ring; additional phosphorylation required"
    
    # --- Step 3: Check for bridging phosphate(s) ---
    # For each phosphate attached to the ring, we check if at least one of its other attached oxygens (not the ring oxygen)
    # is connected to a carbonyl group (an ester linkage, indicating connection with a diacylglycerol fragment).
    bridging_found = False
    
    # Define a helper function to check if an oxygen is part of an ester: it must be connected to a carbon 
    # that has a double bond to an oxygen.
    def is_ester_oxygen(oxygen):
        for nbr in oxygen.GetNeighbors():
            # Skip if neighbor is phosphorus or hydrogen.
            if nbr.GetAtomicNum() in [1,15]:
                continue
            if nbr.GetAtomicNum() == 6:
                # Check bonds of the carbon for a double bond to oxygen.
                for bond in nbr.GetBonds():
                    # Look only for bonds not involving the oxygen we are testing.
                    if oxygen.GetIdx() in [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]:
                        continue
                    if bond.GetBondTypeAsDouble() == 2 and (bond.GetOtherAtom(nbr).GetAtomicNum() == 8):
                        return True
        return False

    for pho_idx in phosphate_set:
        phosphate = mol.GetAtomWithIdx(pho_idx)
        ring_connected_oxygens = set(phosphate_ring_oxygens.get(pho_idx, []))
        # Iterate over all oxygen neighbors of the phosphate.
        for nbr in phosphate.GetNeighbors():
            if nbr.GetAtomicNum() != 8:
                continue
            # Skip if this oxygen is the one that connects to the inositol ring.
            if nbr.GetIdx() in ring_connected_oxygens:
                continue
            # If this oxygen qualifies as part of an ester (i.e. connected to a carbonyl), mark bridging.
            if is_ester_oxygen(nbr):
                bridging_found = True
                break
        if bridging_found:
            break

    if not bridging_found:
        return False, "No bridging phosphate connecting the diacylglycerol fragment to the inositol ring was found"
    
    # If we reach this point, the molecule shows the proper inositol ring, has multiple phosphates (one of which is bridging).
    return True, ("Molecule has an inositol-like ring with sufficient OH-like substituents, multiple phosphate attachments, "
                  "and at least one of the phosphate groups bridges to a diacylglycerol fragment via an ester linkage – "
                  "classified as a phosphoinositide.")

# For quick testing (these lines may be removed or commented out in production)
if __name__ == '__main__':
    # Example SMILES for PIP(18:0/16:0)
    test_smiles = "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O"
    result, reason = is_phosphoinositide(test_smiles)
    print(f"Result: {result}\nReason: {reason}")