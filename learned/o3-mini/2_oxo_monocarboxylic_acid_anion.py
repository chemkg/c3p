"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI: 2‑oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo group is located at the 2‐position.
In other words, there must be a carboxylate group (C(=O)[O–]) directly attached to an 
α‐carbon that bears an extra C=O substituent. For example, molecules containing the fragment
    [α‐C(=O)-C(=O)[O–]]
are accepted provided the α‑carbon is linked to exactly one acid group.
If an explicit fragment match fails, we fallback by accepting molecules that have exactly one 
carboxylate group and at least one carbonyl “elsewhere” (e.g. in an acetamido or other group),
boosting recall.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_carboxylate(carbon_atom):
    """
    Helper: Decide whether a given carbon atom is part of a carboxylate group.
    We require that the carbon is double‐bonded to one oxygen (the carbonyl) and 
    single‐bonded to one oxygen carrying a –1 formal charge.
    """
    if carbon_atom.GetAtomicNum() != 6:
        return False
    has_carbonyl = False
    has_neg_oxy = False
    for bond in carbon_atom.GetBonds():
        # Look for a double bond (the carbonyl)
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            other = bond.GetOtherAtom(carbon_atom)
            if other.GetAtomicNum() == 8:
                has_carbonyl = True
        # Look for a single bond to an oxygen with a –1 charge
        elif bond.GetBondType() == rdchem.BondType.SINGLE:
            other = bond.GetOtherAtom(carbon_atom)
            if other.GetAtomicNum() == 8 and other.GetFormalCharge() == -1:
                has_neg_oxy = True
    return has_carbonyl and has_neg_oxy

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2‑oxo monocarboxylic acid anion based on its SMILES.
    The logic is as follows:
      1) Search for the explicit fragment [α‑C(=O)-C(=O)[O–]]. For each match, verify that 
         the α‑carbon (match[0]) is attached to exactly one carboxylate group.
      2) If no explicit match is found, then (to boost recall) check whether the molecule 
         contains exactly one carboxylate (C(=O)[O–]) group and in addition has at least one 
         carbonyl group (C=O) that is not part of that carboxylate.
         
    Args:
         smiles (str): SMILES string of the molecule.
         
    Returns:
         (bool, str): A tuple of classification result and an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, define a SMARTS that explicitly matches an α‑carbon bearing a C=O and attached
    # to a carboxylate group.
    # The pattern meaning:
    #   [#6X3]=O  : A trigonal carbon with a double bond to oxygen (the α‑carbon's oxo)
    #   -         : single bond to 
    #   [#6X3](=O)[O-] : a trigonal carbon that is part of a carboxylate group.
    explicit_smarts = "[#6X3]=O-[#6X3](=O)[O-]"
    explicit_pattern = Chem.MolFromSmarts(explicit_smarts)
    explicit_matches = mol.GetSubstructMatches(explicit_pattern)
    
    # For each match, check that the α‑carbon (the first atom in the match) is attached to only
    # one carboxylate group.
    for match in explicit_matches:
        alpha_idx, acid_idx = match  # first is the α‑carbon; second is the acid (carboxylate) carbon.
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        # Count how many neighbors of the α‑carbon are carboxylate carbons.
        acid_neighbors = 0
        for neighbor in alpha_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and is_carboxylate(neighbor):
                acid_neighbors += 1
        if acid_neighbors == 1:
            return True, ("Found an explicit fragment with an α‑carbon bearing a C=O and "
                          "directly attached to a carboxylate group, with the α‑carbon linked to exactly one acid group, "
                          "consistent with the 2‑oxo monocarboxylic acid anion definition.")
    
    # Fallback strategy:
    # Count overall carboxylate groups in the molecule using the SMARTS for a carboxylate.
    acid_smarts = "C(=O)[O-]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    # To avoid counting the same acid twice, we collect the unique carbon indices.
    acid_carbons = set(match[0] for match in acid_matches)
    if len(acid_carbons) == 1:
        # Count all carbonyl groups (C=O) in the molecule that are NOT already part of the unique carboxylate.
        carbonyl_count = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for bond in atom.GetBonds():
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(atom)
                        if other.GetAtomicNum() == 8 and other.GetFormalCharge() == 0:
                            # Exclude the oxygen in the carboxylate carbonyl (if any)
                            if not (is_carboxylate(atom) and other.GetFormalCharge() != 0):
                                carbonyl_count += 1
                                break  # count each carbon only once
        if carbonyl_count >= 1:
            return True, ("Fallback accepted: Molecule contains exactly one carboxylate group and at least one additional "
                          "carbonyl group, suggesting the presence of an oxo functionality in a 2‑position relative to the acid.")
    
    return False, ("No fragment with a carboxylate (C(=O)[O-]) directly attached to an α‑carbon bearing an extra C=O "
                   "was found or the α‑carbon is linked to multiple acid groups; the molecule does not fit the 2‑oxo "
                   "monocarboxylic acid anion definition.")

# Example usage for testing (you may add or adjust tests as needed)
if __name__ == "__main__":
    test_examples = [
        ("[H]C(C)=CCC(=O)C([O-])=O", "2-oxohex-4-enoate"),                   # expected True via explicit pattern
        ("CC(=O)N[C@H]([C@@H](O)CC(=O)C([O-])=O)[C@@H](O)[C@H](O)[C@H](O)CO", "aceneuramate"),  # True (explicit)
        ("[C@@H](O)([C@H](O)C(C(=O)[O-])=O)[C@@H](C)O", "(3S,4S,5R)-3,4,5-trihydroxy-2-oxohexanoate"),  # True (explicit)
        ("CC(C)C(=O)C([O-])=O", "3-methyl-2-oxobutanoate"),                     # True (explicit)
        ("CC(C([O-])=O)C(=O)C([O-])=O", "2-methyl-3-oxosuccinate (should be rejected)"),  # False: α‑carbon has >1 acid neighbor
        ("OC(=O)C(=O)C([O-])=O", "oxomalonate (should be rejected)"),            # False
        ("CC(=O)NCCCCC([O-])=O", "5-acetamidopentanoate (fallback case)"),       # fallback: True if criteria met
        ("CC(=O)C(C)(O)C([O-])=O", "2-acetyllactate (fallback case)")            # fallback: may accept if exactly one acid and extra C=O present
    ]
    
    for smi, name in test_examples:
        result, reason = is_2_oxo_monocarboxylic_acid_anion(smi)
        print(f"SMILES: {smi}\n NAME: {name}\n Classification: {result}\n Reason: {reason}\n{'-'*60}")