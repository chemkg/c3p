"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI: 2-oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo group is located at the 2-position.
That is, a carboxylate group (C(=O)[O-]) is directly attached to an alpha-carbon that bears
an additional carbonyl group.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_carboxylate(atom):
    """
    Helper function to decide if a given carbon atom is part of a carboxylate group.
    We require that the carbon (typically sp2) is double-bonded to one oxygen
    (the carbonyl) and single-bonded to another oxygen that carries a -1 formal charge.
    """
    # atom must be a carbon
    if atom.GetAtomicNum() != 6:
        return False
    has_carbonyl = False
    has_oxy_negative = False
    for bond in atom.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() == 8:
                # In a carbonyl group, the oxygen generally is neutral.
                has_carbonyl = True
        elif bond.GetBondType() == rdchem.BondType.SINGLE:
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() == 8 and other.GetFormalCharge() == -1:
                has_oxy_negative = True
    return has_carbonyl and has_oxy_negative

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    It searches for the motif R-C(=O)-C(=O)[O-] in which the carboxylate group is directly attached to
    an alpha-carbon that bears an extra carbonyl (=O). To help avoid false positives (such as in diacid
    moieties where the alpha-carbon connects to more than one acid group), we require that the alpha-carbon
    is attached to exactly one carboxylate group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fits the definition, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS motif for the desired fragment:
    # - [#6] : any carbon (the R-group)
    # - [CX3](=O) : an sp2 carbon bearing an oxo group (the alpha-carbon)
    # - C(=O)[O-] : a carboxylate group
    motif_smarts = "[#6]-[CX3](=O)-C(=O)[O-]"
    motif = Chem.MolFromSmarts(motif_smarts)
    if motif is None:
        return False, "Error in SMARTS pattern"
    
    matches = mol.GetSubstructMatches(motif)
    if not matches:
        return False, ("No fragment of the form R-C(=O)-C(=O)[O-] was found; "
                       "thus no 2-oxo monocarboxylic acid anion motif detected.")
    
    # For each match, verify that the alpha-carbon (match index 1) is directly linked to exactly one carboxylate group.
    # In our SMARTS the atom order is: atom0 (R group), atom1 (alpha-carbon), atom2 (acid carbon of the carboxylate).
    for match in matches:
        # Get the candidate atoms from the match
        r_atom = mol.GetAtomWithIdx(match[0])
        alpha_atom = mol.GetAtomWithIdx(match[1])
        acid_atom = mol.GetAtomWithIdx(match[2])
        
        # Count how many carboxylate carbons are directly attached to the alpha-carbon.
        acid_neighbors = 0
        for neighbor in alpha_atom.GetNeighbors():
            # We are interested in neighbors that are carbon atoms that fulfill the carboxylate criteria.
            if neighbor.GetAtomicNum() == 6 and is_carboxylate(neighbor):
                acid_neighbors += 1
        if acid_neighbors == 1:
            return True, ("Found carboxylate group attached to an alpha carbon that bears a double-bonded "
                          "oxygen (oxo) substituent, with exactly one acid attachment on that alpha carbon, "
                          "fitting the 2-oxo monocarboxylic acid anion definition.")
    
    return False, ("Found a fragment resembling a carboxylate attached to an oxidized alpha-carbon, "
                   "but the alpha carbon is linked to multiple acid groups (or none in the expected sense), "
                   "which does not fit the required 2-oxo monocarboxylic acid anion definition.")

# Example usage for testing (you can add more test cases)
if __name__ == "__main__":
    # Test a couple of examples:
    test_examples = [
        ("[H]C(C)=CCC(=O)C([O-])=O", "2-oxohex-4-enoate"),
        ("CC(=O)N[C@H]([C@@H](O)CC(=O)C([O-])=O)[C@@H](O)[C@H](O)[C@H](O)CO", "aceneuramate"),
        ("CC(C)C(=O)C([O-])=O", "3-methyl-2-oxobutanoate"),
        ("CC(C([O-])=O)C(=O)C([O-])=O", "2-methyl-3-oxosuccinate (false positive expected)"),
        ("CC(=O)NCCCCC([O-])=O", "5-acetamidopentanoate (false negative expected)"),
    ]
    
    for smi, name in test_examples:
        result, reason = is_2_oxo_monocarboxylic_acid_anion(smi)
        print(f"SMILES: {smi}\n NAME: {name}\n Classification: {result}\n Reason: {reason}\n{'-'*60}")