"""
Classifies: CHEBI:18179 phosphoinositide
"""
#!/usr/bin/env python
"""
Classifies: Phosphoinositide
Definition:
  A phosphoinositide is any phosphatidylinositol that is phosphorylated
  at one or more of the hydroxyl groups of the inositol ring.
  
Our revised strategy:
  1. Parse the SMILES and add explicit hydrogens.
  2. Check that the molecule contains at least two ester groups (to capture the diacylglycerol moiety).
  3. Identify a candidate inositol ring:
       • Look for a six‐membered nonaromatic ring where every atom is carbon,
         is sp3 hybridized, and that carries “OH‐like” substituents. We now consider an –OH
         or an oxygen bound to a phosphorus as representing the inositol substituent.
       • We require that at least 4 of the 6 ring carbons have such substituents.
  4. For the candidate inositol ring, count the number of unique phosphorus atoms attached via an oxygen.
       In a phosphatidylinositol (without extra phosphorylation) there will be one bridging phosphate.
       In a phosphoinositide at least one additional phosphate will be present (i.e. count ≥2).
  5. If all conditions are met, we classify as a phosphoinositide.
  
Note: This is a heuristic; there are edge cases and diverse representations.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is defined as a phosphatidylinositol that carries, in addition
    to the phosphate linking group, at least one extra phosphate on the inositol ring.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule is classified as phosphoinositide, False otherwise.
        str: The reasoning behind the classification.
    """
    # Try to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens for correct evaluation of OH groups.
    mol = Chem.AddHs(mol)
    
    # --- Criterion 1: Check existence of at least 2 ester groups ---
    # We use a simple SMARTS pattern to capture an ester bond: carbonyl carbon attached to an oxygen and then a carbon.
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester group(s); at least 2 are required for diacylglycerol chains"

    # --- Criterion 2: Identify an inositol-like ring ---
    # We examine all rings in the molecule.
    ring_found = False
    inositol_ring_atoms = None
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Look for six-membered rings.
        if len(ring) != 6:
            continue
        # Check each atom: must be carbon and sp3 (to avoid aromatic rings)
        valid_ring = True
        substitution_count = 0  # Count of OH-like substituents on ring atoms.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                valid_ring = False
                break
            # Check sp3 hybridization
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                valid_ring = False
                break
            # Look for substituents (neighbors not in the ring) that are “OH-like”
            # An OH substituent is defined either as an oxygen bound to at least one hydrogen,
            # or as an oxygen that in turn is bound to a phosphorus (phosphorylated).
            found_substituent = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen is -OH: either it has a hydrogen neighbor...
                    h_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() == 1]
                    if h_neighbors:
                        found_substituent = True
                        break
                    # ...or if it is bound to a phosphorus (which indicates phosphorylation)
                    p_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() == 15]
                    if p_neighbors:
                        found_substituent = True
                        break
            if found_substituent:
                substitution_count += 1
        # We require at least 4 of the 6 ring carbons have OH-like substitution.
        if valid_ring and substitution_count >= 4:
            ring_found = True
            inositol_ring_atoms = ring
            break

    if not ring_found:
        return False, "No inositol-like cyclohexane ring with sufficient OH-like substituents was found"

    # --- Criterion 3: Check for extra phosphorylation on the inositol ring ---
    # We now iterate over atoms in the putative inositol ring and gather phosphorus atoms attached via oxygen.
    phosphate_neighbors = set()
    for idx in inositol_ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # Only consider neighbors that are O atoms (the bond through which phosphorylation occurs).
            if nbr.GetAtomicNum() != 8:
                continue
            # For this oxygen, check if it is bonded to a phosphorus.
            for subnbr in nbr.GetNeighbors():
                if subnbr.GetAtomicNum() == 15:
                    phosphate_neighbors.add(subnbr.GetIdx())
    # For a plain phosphatidylinositol the inositol ring is linked by one phosphate.
    # A phosphoinositide must have at least one additional phosphate on the ring,
    # hence we require at least 2 unique phosphorus neighbors.
    if len(phosphate_neighbors) < 2:
        return False, f"Only {len(phosphate_neighbors)} phosphate group(s) found on the inositol ring; additional phosphorylation required"
    
    # If all criteria are met, classify as phosphoinositide.
    return True, "Molecule has two ester groups, an inositol-like ring with sufficient substitution, and extra phosphorylation on the ring – classified as a phosphoinositide"

# For quick testing (these lines may be removed or commented out in a production environment)
if __name__ == '__main__':
    # Example SMILES for PIP(18:0/16:0)
    test_smiles = "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O"
    result, reason = is_phosphoinositide(test_smiles)
    print(f"Result: {result}\nReason: {reason}")