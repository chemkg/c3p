"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
Definition: A carbohydrate derivative in which one or more of the oxygens or hydroxy groups 
of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen or any group.
"""

from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string using an improved heuristic.
    
    Heuristic approach:
      1. The molecule must contain at least one sulfur atom.
      2. Search for candidate sugar rings: rings of 5 or 6 atoms that are saturated and non-aromatic.
         For any candidate ring:
           - Count the number of exocyclic hydroxyl (-OH) groups (require at least two).
           - Look for the case where a ring heteroatom (normally oxygen in sugars) is replaced by sulfur,
             OR one of the ring carbons has an external sulfur substituent.
      3. If no suitable ring is found, search for acyclic carbohydrate patterns:
         Look for a contiguous aliphatic chain (3–6 carbons) that is decorated with at least two -OH groups 
         and a sulfur substituent.
    
    Args:
        smiles (str): A SMILES string of the molecule.
    
    Returns:
        (bool, str): True with an explanation if the molecule is classified as a thiosugar;
                     otherwise, False with a reason.
    """
    # Parse SMILES and add explicit hydrogens for accurate detection of -OH groups.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Requirement: at least one sulfur atom must be present.
    if not any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms()):
        return False, "No sulfur atom present, so not a thiosugar"
        
    # Retrieve ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Search through rings for a candidate sugar ring.
    for ring in rings:
        if len(ring) not in [5, 6]:
            continue  # Only consider ring sizes typical for sugars.
        # Exclude rings with any aromatic atoms (e.g., thiophenes).
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Count the number of ring carbons (sp3, non-aromatic) as a rough filter.
        sp3_carbons = [mol.GetAtomWithIdx(idx) for idx in ring 
                       if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(idx).GetIsAromatic()]
        if len(sp3_carbons) < 4:
            continue
        
        # Count the number of hydroxyl (-OH) substituents on ring atoms.
        oh_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # Skip atoms that are part of the ring.
                if nbr.GetSymbol() == "O":
                    # Check if this oxygen qualifies as an -OH (i.e. bonded to at least one hydrogen).
                    if any(nei.GetSymbol() == "H" for nei in nbr.GetNeighbors()):
                        oh_count += 1
        # For a sugar-like structure, require at least two hydroxyl groups.
        if oh_count < 2:
            continue
        
        # CASE 1: The ring heteroatom (normally expected to be O) is replaced by S.
        ring_heteros = [mol.GetAtomWithIdx(idx) for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() in (8, 16)]
        symbols = [atom.GetSymbol() for atom in ring_heteros]
        if "S" in symbols:
            return True, "Thiosugar identified: sugar ring containing a sulfur atom replacing a typical ring oxygen."
        
        # CASE 2: A ring carbon has an external sulfur substituent.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # Only consider carbon atoms.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 16:
                    return True, "Thiosugar identified: sugar ring with an external sulfur substituent replacing a hydroxyl group."
    
    # If no candidate ring was found, attempt to detect acyclic carbohydrate patterns.
    atoms = mol.GetAtoms()
    for atom in atoms:
        # Look for sp3 carbon atoms not in any ring.
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            # Expand the local environment: consider contiguous neighboring sp3 carbons (chain).
            chain = [atom]
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
                    chain.append(nbr)
            if 3 <= len(chain) <= 6:
                chain_oh = 0
                chain_has_S = False
                for at in chain:
                    for nbr in at.GetNeighbors():
                        if nbr.GetIdx() in [a.GetIdx() for a in chain]:
                            continue
                        # Check for hydroxyl group
                        if nbr.GetSymbol() == "O" and any(x.GetSymbol() == "H" for x in nbr.GetNeighbors()):
                            chain_oh += 1
                        if nbr.GetAtomicNum() == 16:
                            chain_has_S = True
                if chain_oh >= 2 and chain_has_S:
                    return True, "Thiosugar identified: acyclic carbohydrate derivative with a sulfur substituent."
    
    return False, "No thiosugar substructure found"

# Example usage:
if __name__ == '__main__':
    # Test one of the known thiosugar examples, e.g. desulfosinigrin.
    test_smiles = "S([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO)C(=NO)CC=C"
    result, reason = is_thiosugar(test_smiles)
    print(result, "->", reason)