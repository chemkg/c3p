"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: Azole
Definition: Any heteroarene with a five‐membered aromatic ring containing at least one nitrogen.
Note: Here we allow for fused ring systems and perform a heuristic check to reject molecules that
appear to be peptides (which sometimes have an azole, for example in histidine).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined as a molecule containing a five‐membered aromatic ring with at least one nitrogen.
    (Fused systems are permitted.) To lower the false-positive rate we apply a peptide heuristic:
    if the molecule has three or more amide bonds and a high molecular weight, we assume it is a peptide
    having an azole side chain.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an azole, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule, which assigns aromaticity, etc.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Error during sanitization: {str(e)}"
    
    # Re-canonicalize the molecule so that RDKit reassigns aromaticity flags properly.
    try:
        can_smiles = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(can_smiles)
        Chem.SanitizeMol(mol)
    except Exception as e:
        # If anything goes wrong, we continue with the original molecule.
        pass

    # Retrieve ring information from the molecule.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    azole_found = False
    # Check through each ring for a five-membered aromatic ring with at least one nitrogen.
    for ring in rings:
        if len(ring) != 5:
            continue
        # Check that every atom in the ring is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Check that at least one atom is nitrogen.
        if not any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
            continue
        # If reached here, we have a candidate five-membered azole ring.
        azole_found = True
        break
    
    if not azole_found:
        return False, "No qualifying five-membered aromatic ring containing nitrogen found"
    
    # Heuristic: try to reject molecules that appear to be peptides.
    # We count amide bonds using a simple SMARTS pattern.
    amide_smarts = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    num_amides = len(amide_matches)
    num_heavy = mol.GetNumHeavyAtoms()
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # If there are three or more amide bonds and the molecular weight is high,
    # then we suspect that the azole ring could be coming from a peptide backbone.
    if num_amides >= 3 and mol_wt > 300:
        return False, "Molecule appears to be a peptide with an azole side chain"
    
    return True, "Found a five-membered aromatic ring containing nitrogen (azole) in the molecule"

# Example usage when running as a script:
if __name__ == "__main__":
    test_smiles = [
        "[N+](C)([C@H](C(=O)[O-])CC=1NC(SC)=NC1)(C)C",  # S-methyl-L-ergothioneine: should be True
        "O(C(=O)C=1NC=CC1)C",                           # Methyl 1H-pyrrole-2-carboxylate: should be True
        "C[Si](Cn1cncn1)(c1ccc(F)cc1)c1ccc(F)cc1",       # flusilazole: should be True
        "C1C=CN=C1",                                   # 3H-pyrrole – note input not in lowercase; canonicalization should help.
        "CCCN(CCOc1c(Cl)cc(Cl)cc1Cl)C(=O)n1ccnc1",       # prochloraz: should be True
        # A likely peptide example (multiple amide bonds and high MW):
        "O=C(N[C@@H](CC=1NC=NC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)[C@H](O)CC2=CC=CC=C2)",
    ]
    
    for s in test_smiles:
        result, reason = is_azole(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")