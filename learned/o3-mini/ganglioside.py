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
      1. At least one sialic acid moiety (e.g., Neu5Ac) as detected by a SMARTS pattern.
      2. A ceramide component characterized by an amide bond and a long aliphatic chain.
      3. An oligosaccharide (sugar) portion as indicated by a ring of 5+ atoms with ≥2 oxygens.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a ganglioside, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check 1: Sialic Acid Moiety ---
    # A simplified SMARTS for a sialic acid (e.g., Neu5Ac) fragment.
    sialic_smarts = "[C@@H]([C@H](O)[C@@H](NC(C)=O))C(=O)O"
    sialic_query = Chem.MolFromSmarts(sialic_smarts)
    if sialic_query is None:
        return False, "Invalid SMARTS pattern for sialic acid"
    if not mol.HasSubstructMatch(sialic_query):
        return False, "No sialic acid moiety detected"
    
    # --- Check 2: Ceramide Component ---
    # Look for an amide bond as a proxy for the ceramide's acyl linkage.
    amide_smarts = "[CX3](=O)[NX3]"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    if amide_query is None:
        return False, "Invalid SMARTS pattern for amide bond (ceramide portion)"
    if not mol.HasSubstructMatch(amide_query):
        return False, "No amide bond found; ceramide (glycosphingolipid) portion missing"
    
    # Also check for a long aliphatic chain.
    # Using a simpler SMARTS: at least 8 consecutive methylene groups.
    chain_smarts = "[CH2]CCCCCCC"  # CH2 followed by 7 carbon atoms
    chain_query = Chem.MolFromSmarts(chain_smarts)
    if chain_query is None:
        return False, "Invalid SMARTS pattern for aliphatic chain"
    if not mol.HasSubstructMatch(chain_query):
        return False, "No long aliphatic chain detected; ceramide portion may be missing"
    
    # --- Check 3: Oligosaccharide Portion ---
    # A crude proxy: at least one ring with 5+ atoms that includes at least 2 oxygen atoms.
    ring_info = mol.GetRingInfo()
    found_sugar = False
    for ring in ring_info.AtomRings():
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if len(ring) >= 5 and oxy_count >= 2:
            found_sugar = True
            break
    if not found_sugar:
        return False, "No sugar ring detected; oligosaccharide portion missing"
    
    # If all three features are found then we classify the molecule as a ganglioside.
    return True, "Molecule contains sialic acid, a ceramide component (amide bond with long chain), and sugar rings; consistent with ganglioside."

# Example usage:
if __name__ == "__main__":
    # Example test from one of the provided cases (you can swap the SMILES string with any from the list)
    test_smiles = "O1[C@](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)C(O)[C@@H](O[C@@H]3CO)OC[C@H](NC(=O)CCCCCCCCCCCCCCC)[C@H](O)/C=C/CCCCCCCCCCCCC)(C[C@H](O)[C@@H](NC(=O)C)C1[C@H](O)[C@H](O)CO)C(O)=O"
    result, reason = is_ganglioside(test_smiles)
    print("Is ganglioside:", result)
    print("Reason:", reason)