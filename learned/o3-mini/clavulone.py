"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone family – esterified prostanoids from marine corals
Definition: A class of esterified prostanoids obtained from marine corals.
Criteria (heuristic):
  • Must contain a cyclopentenone ring (a five-membered ring with a carbonyl group and at least one C=C double bond).
  • Must contain at least one ester group (–O–C(=O)–).
Note: Due to the structural complexity of these natural products, the SMARTS patterns here are only an approximation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone (an esterified prostanoid from marine corals)
    based on its SMILES string.
    
    Checks (heuristically):
    1. The molecule can be parsed.
    2. It contains a cyclopentenone ring (5-membered ring with one ketone and at least one C=C).
    3. It contains at least one ester group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a clavulone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS for a cyclopentenone ring.
    # This pattern looks for a five-membered ring of carbons (R5) where one position is a carbonyl (C=O)
    # and at least one double bond is present.
    cyclopentenone_smarts = "[#6;R5]1[#6;R5]=[#6;R5][#6;R5](=O)[#6;R5]1"
    cyclopentenone_pattern = Chem.MolFromSmarts(cyclopentenone_smarts)
    if not mol.HasSubstructMatch(cyclopentenone_pattern):
        return False, "No cyclopentenone ring found (required prostanoid core missing)"

    # Define a SMARTS for an ester group.
    # This pattern matches an oxygen atom connected to a carbonyl group that in turn is bound to any carbon.
    ester_smarts = "[OX2][CX3](=O)[#6]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester group found (molecule is not esterified)"

    # (Optional) Additional checks may be added. For example, one could check the number of rings,
    # the presence of long alkyl chains, or the molecular weight to further restrict the classification.
    # Here we check that the molecular weight is in a plausible range (>300 Da) for these natural products.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical clavulone"

    return True, "Molecule contains a cyclopentenone ring and at least one ester group consistent with clavulones"

# Example of how one might call the function:
if __name__ == "__main__":
    test_smiles = "BrC1=C[C@@](O)([C@@]([C@@H](OC(=O)C)/C=C/CCCC(OC)=O)(C1=O)[H])C/C=C\\CCCCC"
    result, reason = is_clavulone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)