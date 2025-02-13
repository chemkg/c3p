"""
Classifies: CHEBI:28892 ganglioside
"""
"""
Classifies: Ganglioside
Definition: A molecule composed of a glycosphingolipid (ceramide and oligosaccharide)
with one or more sialic acids linked on the sugar chain.
"""

from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    For our purposes, we require the molecule to contain:
      1. A sialic acid moiety (e.g. Neu5Ac) – approximated using a SMARTS pattern.
      2. A ceramide component – characterized here by first an amide bond 
         (the acyl linkage) and a long aliphatic chain (proxy for fatty acid/sphingosine).
      3. An oligosaccharide (sugar) portion – detected by the presence of sugar rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a ganglioside, False otherwise.
        str: A reason for the classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check 1: Sialic Acid Moiety ---
    # A very simplistic SMARTS for a sialic acid fragment (e.g., Neu5Ac)
    # This pattern looks for a fragment with a carboxylic acid group and an N-acetyl group on an asymmetric center.
    sialic_smarts = "[C@@H]([C@H](O)[C@@H](NC(C)=O))C(=O)O"
    sialic_query = Chem.MolFromSmarts(sialic_smarts)
    if not mol.HasSubstructMatch(sialic_query):
        return False, "No sialic acid moiety detected"
    
    # --- Check 2: Ceramide Component ---
    # In a ceramide there is an amide bond linking the fatty acid and sphingoid bases.
    # We use a SMARTS pattern for an amide "[CX3](=O)[NX3]".
    amide_smarts = "[CX3](=O)[NX3]"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide_query):
        return False, "No amide bond found; ceramide (glycosphingolipid) portion missing"
    
    # Also check for the presence of a long aliphatic chain.
    # We approximate this by requiring a chain of at least 8 consecutive sp3 carbons.
    chain_smarts = "[CX4&H2]{8,}"
    chain_query = Chem.MolFromSmarts(chain_smarts)
    if not mol.HasSubstructMatch(chain_query):
        return False, "No long aliphatic chain detected; ceramide portion may be missing"
    
    # --- Check 3: Oligosaccharide Portion ---
    # Rather than trying to match every type of sugar, we search the molecule’s rings
    # for at least one that is likely to be a sugar. For example, a ring of 5+ atoms that
    # contains at least two oxygens is a crude proxy.
    ring_info = mol.GetRingInfo()
    found_sugar = False
    for ring in ring_info.AtomRings():
        # Count the number of oxygen atoms in the ring
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if len(ring) >= 5 and oxy_count >= 2:
            found_sugar = True
            break
    if not found_sugar:
        return False, "No sugar ring detected; oligosaccharide portion missing"
    
    # If all three features are found then we assume the molecule is a ganglioside.
    return True, "Molecule contains a ceramide part linked to an oligosaccharide with sialic acid moieties; consistent with a ganglioside structure"

# Example usage:
if __name__ == "__main__":
    # (You can replace the SMILES below with one of the provided examples.)
    test_smiles = "O1[C@](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)C(O)[C@@H](O[C@@H]3CO)OC[C@H](NC(=O)CCCCCCCCCCCCCCC)[C@H](O)/C=C/CCCCCCCCCCCCC)(C[C@H](O)[C@@H](NC(=O)C)C1[C@H](O)[C@H](O)CO)C(O)=O"
    result, reason = is_ganglioside(test_smiles)
    print("Is ganglioside:", result)
    print("Reason:", reason)