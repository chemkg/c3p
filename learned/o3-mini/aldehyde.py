"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: Aldehyde (RC(=O)H)
Definition: A compound having a formyl group, that is, a carbonyl (C=O) bonded to at least one hydrogen.
In our improved method we require that the carbonyl carbon is terminal – it is bonded to one heavy atom (the R group) 
or none (in the case of formaldehyde) and at least one hydrogen. This helps reduce false positives.
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An (organic) aldehyde normally has a formyl group where the carbonyl carbon (C) is double bonded to oxygen (O)
    and its other (single) bond is to either one R group or, in the case of formaldehyde, no heavy-atom R group.
    Additionally, at least one hydrogen must be attached to that carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an aldehyde, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Add explicit hydrogens so that we can count them accurately.
    mol = Chem.AddHs(mol)
    
    # Loop over all atoms; we are looking for a carbon candidate:
    # (1) It must be a carbon atom.
    # (2) It must be double bonded to exactly one oxygen (the carbonyl oxygen).
    # (3) Excluding that carbonyl oxygen, the carbon should have at most one heavy atom neighbor (zero for formaldehyde,
    #     one for a typical R-CHO) and no other heavy neighbors that are oxygen (which would indicate acids, esters, etc).
    # (4) It must have at least one hydrogen attached.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        
        # Find double-bonded oxygens.
        carbonyl_os = []
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetSymbol() == "O":
                    carbonyl_os.append(nbr)
        if len(carbonyl_os) != 1:
            # Either no carbonyl or more than one =O (or the carbonyl is ambiguous)
            continue
        
        # Count heavy neighbors other than the carbonyl oxygen.
        heavy_neighbors = []
        for nbr in atom.GetNeighbors():
            # Skip hydrogens and the carbonyl oxygen we already identified.
            if nbr.GetIdx() in [o.GetIdx() for o in carbonyl_os]:
                continue
            if nbr.GetAtomicNum() > 1:  # heavy atom
                heavy_neighbors.append(nbr)
        
        # For a typical aldehyde, the carbonyl carbon should be terminal (only one nonoxygen heavy neighbor).
        # Formaldehyde (H2C=O) would have no such heavy neighbor.
        # Also, if the single neighbor exists, it should not be oxygen.
        if len(heavy_neighbors) == 1 and heavy_neighbors[0].GetSymbol() == "O":
            # This likely indicates a carboxylic acid (or derivative).
            continue
        if not (len(heavy_neighbors) == 0 or len(heavy_neighbors) == 1):
            continue
        
        # Finally, check that at least one hydrogen is attached.
        # Using GetTotalNumHs(), which counts both explicit and implicit hydrogens.
        if atom.GetTotalNumHs() < 1:
            continue
        
        # If we have reached here then we consider that a valid aldehyde group is present.
        return True, ("Aldehyde group detected: carbonyl (C=O) "
                      "bonded to at least one hydrogen and a terminal substituent "
                      "(R group or none, as in formaldehyde).")
    
    return False, "No aldehyde functional group (RC(=O)H) detected."

# Example usage (uncomment to test):
# test_smiles = "OC=1C(=CC=CC1)/C=C/C=O"  # (E)-3-(2-Hydroxyphenyl)-2-propenal
# result, reason = is_aldehyde(test_smiles)
# print(result, reason)