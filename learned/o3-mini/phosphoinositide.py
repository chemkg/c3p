"""
Classifies: CHEBI:18179 phosphoinositide
"""
#!/usr/bin/env python
"""
Classifies: Phosphoinositide

Definition:
  A phosphoinositide is any phosphatidylinositol that is phosphorylated
  at one or more of the hydroxyl groups of the inositol ring.
  
Revised strategy:
  1. Parse the SMILES and add explicit hydrogens.
  2. Verify that the molecule contains two or more ester groups
     (to capture the diacylglycerol fragment). Also record the ester oxygen atoms.
  3. Identify a candidate inositol-like ring:
       • Look for a six-membered ring in which every atom is a carbon and is sp3 hybridized.
       • For each ring carbon, check for at least one substituent on an external oxygen that is “OH‐like.”
         An “OH‐like” substituent is defined as an oxygen (with formal charge 0) that has at least one hydrogen,
         or an oxygen bound to phosphorus.
       • Accept the ring if at least 4 of its 6 carbons carry an OH-like substituent.
  4. For the candidate inositol ring, look over its substituents. Gather all phosphate (P) atoms that are attached
     via an oxygen. In addition, we require that at least one phosphate be “bridging” to the diacylglycerol,
     as indicated by its oxygen being one of the ester oxygens found earlier.
  5. Finally, to be a phosphoinositide we require that the ring shows at least two unique phosphate attachments
     (one for the glycerol bridge and at least one extra phosphorylation).
     
Note:
  This approach is heuristic; depending on non‐standard representations some edge cases may be missed.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is defined as a phosphatidylinositol that (in addition to its bridging phosphate)
    carries at least one extra phosphate on the inositol ring.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phosphoinositide, False otherwise.
        str: The reasoning behind the classification.
    """
    # Attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that OH groups are correctly represented.
    mol = Chem.AddHs(mol)
    
    # --- Criterion 1: Check for the presence of diacylglycerol ester groups ---
    # Use a simple SMARTS to find an ester motif: carbonyl carbon attached to an oxygen which is then attached to a carbon.
    ester_smarts = "C(=O)O[C]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester group(s); at least 2 are required to indicate diacylglycerol chains"
    
    # Record the indices of ester oxygen atoms (the atom at position 1 in the SMARTS pattern).
    ester_oxygens = set(match[1] for match in ester_matches)
    
    # --- Criterion 2: Identify an inositol-like ring.
    # We iterate over all rings; candidate rings must be six-membered, all carbon, sp3.
    ring_info = mol.GetRingInfo()
    candidate_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        valid_ring = True
        substitution_count = 0  # count of ring carbons with an OH-like substituent
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Must be carbon and sp3 hybridized (to avoid aromatic or hetero rings)
            if atom.GetAtomicNum() != 6 or atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                valid_ring = False
                break
            # Look among neighbors not in the ring. We consider a substituent valid if it is
            # an oxygen (with formal charge 0) that carries at least one hydrogen OR is bound to phosphorus.
            found_substituent = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check for an OH group: must have at least one hydrogen and formal charge zero.
                    if nbr.GetFormalCharge() == 0 and any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                        found_substituent = True
                        break
                    # Or, check if the oxygen is bound to phosphorus (even if deprotonated).
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
    
    # --- Criterion 3: Count phosphate attachments on the candidate inositol ring ---
    # For each ring carbon, look at neighbor oxygen atoms and check for attached phosphorus.
    phosphate_indices = set()
    bridging_found = False  # A flag to confirm the phosphate linking the ring to a diacylglycerol (i.e. via an ester oxygen)
    for idx in candidate_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # Skip if neighbor is part of the ring.
            if nbr.GetIdx() in candidate_ring:
                continue
            if nbr.GetAtomicNum() != 8:
                continue
            # For each oxygen substituent, see if it is bonded to a phosphorus.
            for subnbr in nbr.GetNeighbors():
                if subnbr.GetAtomicNum() == 15:
                    phosphate_indices.add(subnbr.GetIdx())
                    # Check if this oxygen is one of the ester oxygens (bridging the diacylglycerol).
                    if nbr.GetIdx() in ester_oxygens:
                        bridging_found = True
    # A normal phosphatidylinositol will always have one bridging phosphate;
    # a phosphoinositide should have at least one extra phosphate on the ring (i.e. count at least 2).
    if len(phosphate_indices) < 2:
        return False, f"Only {len(phosphate_indices)} phosphate group(s) found on the inositol ring; additional phosphorylation required"
    if not bridging_found:
        return False, "No bridging phosphate connecting the diacylglycerol fragment to the inositol ring was found"
    
    return True, ("Molecule has diacylglycerol ester groups, an inositol-like ring with sufficient OH-like substituents, "
                  "and extra phosphorylation on the ring (with a bridging phosphate) – classified as a phosphoinositide")

# For quick testing (these lines may be removed or commented out in production)
if __name__ == '__main__':
    # Example SMILES for PIP(18:0/16:0)
    test_smiles = "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O"
    result, reason = is_phosphoinositide(test_smiles)
    print(f"Result: {result}\nReason: {reason}")