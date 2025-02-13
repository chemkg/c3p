"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: secondary alcohol
Definition: A secondary alcohol is a compound in which a hydroxy group (-OH)
is attached to a saturated carbon atom that has exactly two carbon substituents.
This function uses RDKit to parse the SMILES string, add explicit hydrogens, and 
then looks for an oxygen (from an -OH group) attached to a sp3 carbon meeting the criteria.
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule contains at least one secondary alcohol group.
    
    A secondary alcohol is defined as a hydroxyl group (-OH) attached to an sp3 (tetrahedral) 
    carbon atom that has exactly two other carbon substituents and one hydrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if at least one secondary alcohol group is found, False otherwise.
        str: Explanation or reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for accurate count.
    mol = Chem.AddHs(mol)
    
    # Loop over atoms looking for oxygen atoms possibly from an alcohol group.
    for atom in mol.GetAtoms():
        # Check if the atom is oxygen.
        if atom.GetAtomicNum() != 8:
            continue
        neighbors = atom.GetNeighbors()
        # For an -OH group the oxygen should be connected to exactly one hydrogen and one carbon.
        if len(neighbors) != 2:
            continue

        # Identify which neighbor is hydrogen and which is carbon.
        oxygen_neighbor_H = None
        oxygen_neighbor_C = None
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 1:
                oxygen_neighbor_H = nbr
            elif nbr.GetAtomicNum() == 6:
                oxygen_neighbor_C = nbr
        if oxygen_neighbor_H is None or oxygen_neighbor_C is None:
            continue  # Not an -OH group.

        # Now, examine the carbon attached to this oxygen.
        carbon_atom = oxygen_neighbor_C
        # Ensure the carbon is sp3 (tetrahedral).
        if carbon_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue

        # Count the number of carbon neighbors attached to the carbon atom, excluding the oxygen.
        carbon_neighbors = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 2:
            continue

        # Count hydrogen atoms attached to the carbon.
        hydrogen_count = sum(1 for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
        if hydrogen_count != 1:
            continue

        # If all conditions are satisfied, we've found a secondary alcohol group.
        return True, f"Molecule contains a secondary alcohol group at carbon atom index {carbon_atom.GetIdx()}."
    
    # No matching secondary alcohol group was found.
    return False, "No secondary alcohol group found."

# Example usage (remove or modify these lines for integration in your pipeline):
if __name__ == "__main__":
    # A short list of examples (from the provided list) for testing:
    examples = [
        ("O[C@H]([C@H](O)CO)CCO", "2-deoxy-D-erythro-pentitol"),
        ("C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@](C)(CCC=C(C)C)O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@]3(C)CC[C@H](O)C(C)(C)[C@H]23)[C@H](O)[C@H](O)[C@H]1O", "ginsenoside Re"),
        ("O[C@@H](CC\\C=C\\c1ccccc1)CCc1ccc(O)c(O)c1", "(-)-(3S)-1-(3,4-dihydroxyphenyl)-7-phenyl-(6E)-6-hepten-3-ol"),
        ("C[C@H](O)[C@H](C)O", "(S,S)-butane-2,3-diol")
    ]
    for smi, name in examples:
        flag, reason = is_secondary_alcohol(smi)
        print(f"SMILES: {smi}\nName: {name}\nClassification: {flag}\nReason: {reason}\n")