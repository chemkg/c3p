"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A compound that contains exactly two ketone functionalities (diketone)

A diketone is defined here as a molecule in which there exist exactly 2 carbonyl groups 
(i.e. C=O) where the carbonyl carbon is bonded to two carbon atoms (thus excluding aldehydes
or carbonyls bonded to heteroatoms). In addition, we verify that the only carbonyls present 
in the molecule are these two ketones.
"""

from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is defined as a compound containing exactly two ketone groups.
    Here a ketone group is defined as a carbonyl (C=O) where the carbon is bound 
    to two carbons (i.e. excluding aldehydes or carbonyls bonded to heteroatoms).
    
    The approach first parses the SMILES into a molecule and then iterates over its atoms.
    For each carbon we check if it has a double bond to oxygen and has two carbon neighbors
    (aside from that oxygen). We then also count each oxygen in such a C=O bond.
    Finally, we require that no additional C=O groups exist in the molecule; i.e. the total 
    number of carbonyl groups equals the number of ketone groups (which must equal 2).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as a diketone, False otherwise.
        str: Explanation for the result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # sets to hold indices of ketone carbon atoms and their paired oxygen atoms
    ketone_carbons = set()
    ketone_oxygens = set()
    
    # loop over atoms of the molecule
    for atom in mol.GetAtoms():
        # We only examine carbon atoms (atomic number 6)
        if atom.GetAtomicNum() != 6:
            continue
        # get neighbors and check for a double bond to oxygen
        neighbors = atom.GetNeighbors()
        has_double_bound_oxygen = False
        oxygen_idx = None
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 8:
                # check if the bond is double
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    has_double_bound_oxygen = True
                    oxygen_idx = nbr.GetIdx()
                    # We assume only one C=O is considered for a given carbon
                    break
        if not has_double_bound_oxygen:
            continue
        
        # Now check that (aside from the oxygen) the carbon has exactly 2 other carbon neighbors.
        # Count only bonds that are not the C=O bond.
        carbon_neighbor_count = 0
        for nbr in neighbors:
            # skip the oxygen that participates in the double bond
            if oxygen_idx is not None and nbr.GetIdx() == oxygen_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                carbon_neighbor_count += 1
        if carbon_neighbor_count == 2:
            ketone_carbons.add(atom.GetIdx())
            ketone_oxygens.add(oxygen_idx)
    
    # Now count all carbonyl groups in the molecule
    # We define a carbonyl group as any carbon atom double bonded to an oxygen.
    total_carbonyl_carbons = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    total_carbonyl_carbons.add(atom.GetIdx())
                    break

    n_ketones = len(ketone_carbons)
    n_total_carbonyls = len(total_carbonyl_carbons)
    
    # Check the diketone criteria:
    # 1) exactly 2 ketone groups should be found by our definition
    # 2) there should be no extra carbonyl groups in the molecule.
    if n_ketones == 2 and n_total_carbonyls == 2 and len(ketone_oxygens) == 2:
        return True, "Compound contains exactly 2 ketone groups."
    else:
        # Build a detailed message
        if n_ketones != 2:
            return False, f"Compound contains {n_ketones} ketone group(s) (by strict ketone criteria), which does not equal 2."
        elif n_total_carbonyls != 2:
            return False, f"Compound contains {n_total_carbonyls} carbonyl group(s); extra carbonyl functionality detected beyond the 2 ketone groups."
        else:
            return False, "Compound does not meet diketone criteria."

# Example usage (you can remove these lines if integrating into another codebase):
if __name__ == "__main__":
    test_smiles = [
        "O=C(CCCCCCCCCCCCCCC)CC(=O)CCCCCCCCC",  # 10,12-Heptacosanedione (diketone, expected True)
        "CCC(=O)C(C)=O",                        # pentane-2,3-dione (diketone, expected True)
        "O=C(CC(=O)c1ccccc1)c1ccccc1",            # dibenzoylmethane (diketone, expected True)
        "OC1C(C(=O)C(C1=O)C(=O)C(CC)C)CC=C(C)C",  # Adhumulinic acid (3 ketones, expected False)
    ]
    for smi in test_smiles:
        result, reason = is_diketone(smi)
        print(f"SMILES: {smi}")
        print(f"Result: {result}, Reason: {reason}")
        print("-"*50)