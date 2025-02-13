"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI: 2-oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo group is located at the 2‐position.
In other words, there is a carboxylate group (C(=O)[O-]) directly attached to an α‐carbon that 
bears an extra C=O substituent. For example, molecules containing the fragment
    [*]-C(=O)-C(=O)[O-]
(such as glyoxylate, 2-oxohex-4-enoate, aceneuramate, etc.) are accepted.
We further require that the α‑carbon (the one bearing the ketone group) is attached to exactly one 
carboxylate group – to avoid diacid (or polyacid) cases.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_carboxylate(atom):
    """
    Helper: Given a carbon atom, determine if it is part of a carboxylate group.
    We require that the carbon is double‐bonded to one oxygen (carbonyl) and 
    single‐bonded to an oxygen carrying a –1 formal charge.
    """
    if atom.GetAtomicNum() != 6:
        return False
    has_carbonyl = False
    has_neg_oxy = False
    for bond in atom.GetBonds():
        # Check for double bonds (carbonyl)
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() == 8:
                has_carbonyl = True
        # Check for single-bonded oxygens with a -1 formal charge
        elif bond.GetBondType() == rdchem.BondType.SINGLE:
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() == 8 and other.GetFormalCharge() == -1:
                has_neg_oxy = True
    return has_carbonyl and has_neg_oxy

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES.
    The algorithm first finds candidate carboxylate groups (C(=O)[O-]). For each such group, 
    we locate the carbon it is attached to (the candidate alpha‐carbon). If that carbon has 
    (a) a double‐bonded oxygen (that is not the carboxylate group) and 
    (b) is linked to exactly one carboxylate group, then we decide the fragment R–C(=O)–C(=O)[O–]
    is present.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule fits the definition, False otherwise.
       str: Explanation for the decision.
       
    In borderline cases (or if structural features are ambiguous) one could return (None, None).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS to find carboxylate groups explicitly.
    acid_smarts = "C(=O)[O-]"
    acid_mol = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_mol)
    
    if not acid_matches:
        return False, "No carboxylate (C(=O)[O-]) group found"
    
    # For each found carboxylate (we take the carbon of the acid group),
    # check its single (non-oxygen) neighbor, the potential α-carbon.
    for match in acid_matches:
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Get all atoms directly bonded to acid_carbon
        neighbors = acid_carbon.GetNeighbors()
        # We want a neighbor that is a carbon or any atom (including H via implicit valence)
        # In many correct examples (e.g. glyoxylate, 2-oxohex-4-enoate) this join is direct.
        candidate_alphas = [n for n in neighbors if n.GetAtomicNum() == 6]
        if not candidate_alphas:
            continue  # nothing to check for this acid group
        
        # For each candidate alpha, check that it actually bears a C=O (ketone or aldehyde) bond.
        for alpha in candidate_alphas:
            # Count how many carboxylate groups are attached to alpha.
            acid_neighbors = 0
            for nb in alpha.GetNeighbors():
                if nb.GetAtomicNum() == 6 and is_carboxylate(nb):
                    acid_neighbors += 1
            # We require exactly one carboxylate attachment (the one we are considering)
            if acid_neighbors != 1:
                continue

            # Now check that α-carbon bears a double-bonded oxygen (other than the one in the carboxylate).
            has_oxo = False
            for bond in alpha.GetBonds():
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(alpha)
                    # Exclude the acid carbon that is already counted if it is oxygen (it isn’t)
                    if other.GetAtomicNum() == 8:
                        # We allow both ketone and aldehyde oxygen (formal charge usually 0)
                        # This is our “oxo” substituent.
                        has_oxo = True
                        break
            if has_oxo:
                return True, ("Found a carboxylate group whose attached (α-) carbon "
                              "carries a C=O substituent and is linked to exactly one acid group, "
                              "consistent with the 2‑oxo monocarboxylic acid anion definition.")
    # If no candidate passes, explain.
    return False, ("No fragment with a carboxylate (C(=O)[O-]) directly attached to an α‑carbon "
                   "bearing an extra C=O was found; the molecule does not fit the 2‑oxo monocarboxylic acid "
                   "anion definition.")

# Example usage for testing. (You can add your own test cases.)
if __name__ == "__main__":
    test_examples = [
        ("[H]C(C)=CCC(=O)C([O-])=O", "2-oxohex-4-enoate"),  # expected True
        ("CC(=O)N[C@H]([C@@H](O)CC(=O)C([O-])=O)[C@@H](O)[C@H](O)[C@H](O)CO", "aceneuramate"),  # True
        ("[C@@H](O)([C@H](O)C(C(=O)[O-])=O)[C@@H](C)O", "(3S,4S,5R)-3,4,5-trihydroxy-2-oxohexanoate"),  # True
        ("CC(C)C(=O)C([O-])=O", "3-methyl-2-oxobutanoate"),  # True
        ("CC(C([O-])=O)C(=O)C([O-])=O", "2-methyl-3-oxosuccinate (should be rejected)"),  # False: α‐carbon has >1 acid
        ("OC(=O)C(=O)C([O-])=O", "oxomalonate (should be rejected)"),  # False
        ("[H]C(=O)C([O-])=O", "glyoxylate"),  # True: even if the R-group is H
        ("CC(=O)NCCCCC([O-])=O", "5-acetamidopentanoate"),  # May be borderline; here our criteria require a contiguous oxo…
    ]
    
    for smi, name in test_examples:
        result, reason = is_2_oxo_monocarboxylic_acid_anion(smi)
        print(f"SMILES: {smi}\n NAME: {name}\n Classification: {result}\n Reason: {reason}\n{'-'*60}")