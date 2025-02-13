"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: Any nucleoside where the sugar component is D-ribose.
A ribonucleoside is a nucleoside whose sugar is D-ribose – that is, a furanose ring (5-membered) 
with four carbons and one oxygen; the sugar is typically decorated with –OH groups (plus a CH2OH on C5'), 
and the anomeric carbon (C1’) is linked to a nucleobase (an aromatic, N-containing heterocycle).
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside (i.e. a nucleoside where the sugar is D-ribose) 
    based on its SMILES string.
    
    The function uses two SMARTS patterns to try to pick up the ribose moiety. It then checks 
    that one of the ribose ring atoms – the anomeric carbon – is connected to a nitrogen that belongs 
    to an aromatic heterocycle (i.e. the nucleobase).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a ribonucleoside, otherwise False.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
 
    # Define SMARTS patterns for a D-ribose ring as seen in many nucleosides.
    # Note: These patterns assume proper stereochemical descriptors.
    # Pattern 1 (as seen in e.g. adenosine, ribothymidine):
    ribose_smarts1 = "[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"
    # Pattern 2 (sometimes the CH2OH can be drawn as a separate fragment, e.g. dihydrouridine):
    ribose_smarts2 = "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)"
    
    # Try each ribose SMARTS pattern to find a ribose-like motif.
    ribose_matches = []
    for smarts in (ribose_smarts1, ribose_smarts2):
        query = Chem.MolFromSmarts(smarts)
        if query is None:
            continue  # skip if SMARTS pattern cannot be parsed
        matches = mol.GetSubstructMatches(query)
        if matches:
            # For this implementation, we take the first match
            ribose_matches = list(matches)[0]  # a tuple of atom indices for the matching ring
            break
    
    # If no ribose pattern was found:
    if not ribose_matches:
        return False, "No D-ribose sugar moiety detected"
    
    # At this point, we have found a ribose-like ring.
    # Next we want to check for attachment of a nucleobase.
    # The glycosidic (anomeric) carbon in many nucleosides appears as the first atom of the SMARTS match.
    # (For example, in adenosine the anomeric carbon is written as "[C@@H]" which is attached to a nitrogen of the base.)
    anomeric_atom_idx = ribose_matches[0]
    anomeric_atom = mol.GetAtomWithIdx(anomeric_atom_idx)
    
    # Look for a neighbor of the anomeric carbon that is outside the matched ribose ring.
    # We require that the neighbor must be a nitrogen and be aromatic – a common indicator
    # of a nucleobase.
    ribose_atom_set = set(ribose_matches)
    nucleobase_attached = False
    for neighbor in anomeric_atom.GetNeighbors():
        if neighbor.GetIdx() not in ribose_atom_set:
            if neighbor.GetAtomicNum() == 7 and neighbor.GetIsAromatic():
                nucleobase_attached = True
                break

    if not nucleobase_attached:
        return False, "Ribose detected but no attached nucleobase (nitrogen in aromatic ring) found"
    
    return True, "Contains D-ribose sugar moiety with an attached nucleobase"

# For testing (uncomment the following lines):
# test_smiles = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"  # adenosine
# print(is_ribonucleoside(test_smiles))