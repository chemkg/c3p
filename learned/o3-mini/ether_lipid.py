"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: ether lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more of the carbon atoms on glycerol 
             is bonded to an alkyl chain via an ether linkage instead of an ester linkage.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    It first checks if a glycerol backbone is present (using a simple SMARTS pattern) and then verifies
    that at least one alkyl chain is attached via an ether linkage rather than an ester linkage.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an ether lipid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Optional: Sanitize molecule
    Chem.SanitizeMol(mol)
    
    # Check for a glycerol backbone.
    # Here we attempt a very simple SMARTS pattern that matches a glycerol-like fragment.
    # This pattern looks for a 3-carbon chain with hydroxyl (or ether) substituents: OCC(O)CO.
    # Note that many lipid variants replace one or more -OH with substituents but we use this as a minimal check.
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone detected (using pattern OCC(O)CO)"
    
    # Function to decide if an oxygen is part of an ester group.
    def is_ester_oxygen(o_atom):
        # In an ester, the oxygen is connected to a carbon that is double-bonded to an O.
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                # Check if this carbon has a double bond to oxygen
                for bond in nbr.GetBonds():
                    if bond.GetBondTypeAsDouble() == 2.0:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # Now scan molecule for ether oxygen atoms.
    # Ether oxygen in a lipid is an oxygen connected to two carbons (and not to hydrogen) and not involved in a carbonyl (ester linkage).
    ether_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            # Count heavy-atom neighbors (exclude hydrogens)
            heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
            if len(heavy_neighbors) == 2:
                # Skip if this oxygen is an -OH (has hydrogen) by checking explicit hydrogens count.
                # Also, if the oxygen is part of an ester, skip.
                if atom.GetTotalNumHs() > 0:
                    continue
                if is_ester_oxygen(atom):
                    continue
                # Check both neighbors are carbons
                if all(nbr.GetAtomicNum() == 6 for nbr in heavy_neighbors):
                    ether_found = True
                    break
    
    if not ether_found:
        return False, "No non-ester (ether) linkage found attached to glycerol"
    
    return True, "Molecule contains glycerol backbone with at least one alkyl chain attached via an ether linkage"

# Example usage:
if __name__ == "__main__":
    # Test one of the given examples:
    test_smiles = "CCCCCCCCCCCCCC\\C=C\\OC[C@@H](O)CO"  # 1-[(E)-hexadecen-1-yl]-sn-glycerol
    result, reason = is_ether_lipid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)