"""
Classifies: CHEBI:35785 sphingoid
"""
#!/usr/bin/env python
"""
Classifies: Sphingoid (Sphinganine, its homologs and stereoisomers, 
and the hydroxy and unsaturated derivatives of these compounds.)
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule belongs to the sphingoid class.
    Sphingoid compounds (e.g. sphinganine, sphingosine, phytosphingosine and their derivatives)
    are defined by a long aliphatic backbone (typically 14–26 carbons) that includes an
    amino alcohol (or its substituted/unsaturated/deoxy derivative) portion.
    
    This function performs a series of heuristic checks:
      1. It verifies that the molecule has a long (non‐cyclic) alkyl chain.
      2. It requires the presence of at least one nitrogen atom (to allow an amino or acylamino group).
      3. It requires at least one non‐carbonyl –OH (or a deoxy variant).
      4. It searches for a “core” sphingoid fragment.
         For instance, sphinganine has the fragment –[C](O)-[C](N)CO (ignoring stereochemistry).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is likely a sphingoid, False otherwise.
        str: Reason for classification.
    """
    # parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # step 1: check for a long aliphatic (non‐cyclic) chain.
    # We count the number of carbon atoms that are not in any ring.
    non_ring_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if len(non_ring_carbons) < 10:
        return False, "Not enough non‐ring carbon atoms to represent a long aliphatic chain"
    
    # step 2: check for nitrogen presence (note: sphingoid bases have an amino function, free or acylated)
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(nitrogens) == 0:
        return False, "No nitrogen found (expect an amino or acylamino group in sphingoid compounds)"
    
    # step 3: check for at least one hydroxyl (OH) group.
    # We look for atoms “O” with one hydrogen attached (i.e. –OH)
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(oh_pattern):
        # Allow a possibility of a deoxy sphingoid derivative
        deoxy_note = " (no OH found, so may be a deoxy derivative)"
    else:
        deoxy_note = ""
    
    # step 4: search for a sphingoid core motif.
    # Our heuristic is to find two connected carbons where one bears an -OH and the neighboring carbon
    # bears a nitrogen substituent and (in many cases) is attached to a hydroxymethyl (-CO) group.
    # (We ignore stereochemistry.)
    #
    # We define a couple of alternative SMARTS for variations:
    #
    # Pattern A: “normal” sphingoid base: –[C](O)-[C](N)CO
    pattern_A = Chem.MolFromSmarts("[C;!R](O)-[C;!R]([#7])(CO)")
    # Pattern B: a variant where the first carbon may be deoxy: –[C]-[C](N)CO
    pattern_B = Chem.MolFromSmarts("[C;!R]-[C;!R]([#7])(CO)")
    # Pattern C: a variant where the terminal hydroxymethyl is not present: –[C](O)-[C](N)
    pattern_C = Chem.MolFromSmarts("[C;!R](O)-[C;!R]([#7])")
    
    matches_A = mol.GetSubstructMatches(pattern_A)
    matches_B = mol.GetSubstructMatches(pattern_B)
    matches_C = mol.GetSubstructMatches(pattern_C)
    
    if matches_A:
        return True, "Matches normal sphingoid backbone motif" + deoxy_note
    elif matches_B:
        return True, "Matches deoxy/unsaturated sphingoid backbone motif (lacking OH at first carbon)" + deoxy_note
    elif matches_C:
        return True, "Matches partial sphingoid backbone motif (may be a deoxy derivative)" + deoxy_note
    else:
        return False, "No sphingoid backbone motif found in molecule"

# Example usage:
if __name__ == "__main__":
    # Some example SMILES strings from the provided list.
    examples = [
        "CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C",  # glycosylated sphingoid derivative
        "CCCCCCCCCCCC[C@@H](O)[C@@H](O)[C@@H](N)CO",  # C16 phytosphingosine
        "CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H]([NH3+])CO",  # tetradecasphingosine (with unsaturation)
    ]
    
    for smi in examples:
        flag, reason = is_sphingoid(smi)
        print(f"SMILES: {smi}\n  --> Sphingoid: {flag}, Reason: {reason}\n")