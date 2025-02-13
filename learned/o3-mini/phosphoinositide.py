"""
Classifies: CHEBI:18179 phosphoinositide
"""
#!/usr/bin/env python
"""
Classifies: Phosphoinositide 
Definition:
  Any phosphatidylinositol that is phosphorylated at one or more of the hydroxyl groups of inositol.
  
Our strategy:
  1. Parse the SMILES.
  2. Count phosphorus atoms – require at least 2.
  3. Look for a 6–membered (cyclohexane) ring in which at least 4 carbons have an –OH substituent.
     (This is our approximate inositol ring candidate.)
  4. Identify at least 2 ester groups (C(=O)O–) that suggest the diacylglycerol moiety.
  5. If all conditions are met, we classify the molecule as a phosphoinositide.
  
Note: This is an approximate filter; for a very challenging class there may be edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is defined as any phosphatidylinositol that is phosphorylated
    on one or more of the hydroxyl groups of the inositol ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a phosphoinositide, False otherwise
        str: Reason for the classification result
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure hydrogens are computed so that we can check for -OH groups correctly
    mol = Chem.AddHs(mol)

    # 1. Count phosphorus atoms (atomic number 15). We require at least 2:
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 2:
        return False, f"Found only {p_count} phosphorus atom(s); at least 2 required for a phosphoinositide"

    # 2. Look for a candidate inositol ring.
    # We search for a 6-membered nonaromatic ring composed entirely of carbon atoms.
    # Then, we count how many exocyclic -OH groups are attached to these ring carbons.
    ring_found = False
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is carbon
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            continue
        # Count hydroxyl (–OH) substituents on ring atoms.
        oh_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors that are not in the ring
            for nbr in atom.GetNeighbors():
                # We want -OH groups: an oxygen atom bonded to at least one hydrogen.
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                    # Get the total number of explicit hydrogens on this O atom.
                    # (After AddHs, explicit hydrogens are present.)
                    h_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() == 1]
                    if len(h_neighbors) >= 1:
                        oh_count += 1
        # Require that the ring has at least 4 hydroxyl groups (an approximate criterion)
        if oh_count >= 4:
            ring_found = True
            break

    if not ring_found:
        return False, "No inositol-like cyclohexane ring with sufficient hydroxyl groups was found"

    # 3. Look for at least two ester groups. Ester groups are typically described by a carbonyl attached to an oxygen bonded to carbon.
    # We use a SMARTS pattern for an ester: [CX3](=O)O[C;!$(O=C)]
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester group(s); at least 2 required to indicate diacylglycerol chains"

    return True, "Molecule has ≥2 phosphorus atoms, an inositol-like ring, and two esters – classified as a phosphoinositide"

# For quick testing (remove or comment these lines in production)
if __name__ == '__main__':
    # example SMILES for one of the provided phosphoinositide examples:
    test_smiles = "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O"
    result, reason = is_phosphoinositide(test_smiles)
    print(f"Result: {result}\nReason: {reason}")