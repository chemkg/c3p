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
    For our purposes we require the molecule to contain:
      1. At least one sialic acid moiety (e.g., Neu5Ac or Neu5Gc). Rather than a strict SMARTS with chirality,
         we now use two alternative SMARTS patterns that ignore chirality.
      2. A ceramide component characterized by an amide bond and a long aliphatic chain.
      3. An oligosaccharide portion, as indicated by the presence of at least one sugar ring
         (a ring of at least 5 atoms with at least 2 oxygens).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a ganglioside, False otherwise.
        str: Explanation for the classification decision.
    """

    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ---- Check 1: Sialic Acid Moiety ----
    # We use two alternative SMARTS patterns for sialic acid.
    # Note: We remove chirality markers to allow a more flexible match.
    # Pattern 1: Typical Neu5Ac moiety
    sialic_smarts1 = "C1(C(=O)O)C(O)C(O)C(NC(C)=O)C(O)1O"
    # Pattern 2: Neu5Gc variant (note the extra oxygen in the N-acyl group)
    sialic_smarts2 = "C1(C(=O)O)C(O)C(O)C(NC(CO)=O)C(O)1O"
    sialic_query1 = Chem.MolFromSmarts(sialic_smarts1)
    sialic_query2 = Chem.MolFromSmarts(sialic_smarts2)
    # If neither sialic acid pattern is found, then reject the molecule.
    if not (mol.HasSubstructMatch(sialic_query1) or mol.HasSubstructMatch(sialic_query2)):
        return False, "No sialic acid moiety detected"

    # ---- Check 2: Ceramide Component ----
    # We detect an amide bond as part of the ceramide’s acyl linkage.
    amide_smarts = "[CX3](=O)[NX3]"  # general amide group
    amide_query = Chem.MolFromSmarts(amide_smarts)
    if amide_query is None:
        return False, "Invalid SMARTS pattern for amide bond (ceramide portion)"
    if not mol.HasSubstructMatch(amide_query):
        return False, "No amide bond found; ceramide (glycosphingolipid) portion missing"

    # Check for a long aliphatic chain (as found in the fatty acid of the ceramide)
    # Here we require a chain of at least 9 contiguous carbons.
    # The pattern "[CH2]CCCCCCCC" represents one CH2 group plus eight additional carbon atoms.
    chain_smarts = "[CH2]CCCCCCCC"
    chain_query = Chem.MolFromSmarts(chain_smarts)
    if chain_query is None:
        return False, "Invalid SMARTS pattern for long aliphatic chain"
    if not mol.HasSubstructMatch(chain_query):
        return False, "No long aliphatic chain detected; ceramide portion may be missing"

    # ---- Check 3: Oligosaccharide Portion ----
    # As a proxy for sugars we check if the molecule has at least one ring of at least 5 atoms
    # with at least two oxygens.
    ring_info = mol.GetRingInfo()
    found_sugar = False
    for ring in ring_info.AtomRings():
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if len(ring) >= 5 and oxy_count >= 2:
            found_sugar = True
            break
    if not found_sugar:
        return False, "No sugar ring detected; oligosaccharide (sugar) portion missing"

    # All three key features are detected.
    return True, ("Molecule contains at least one sialic acid moiety, an amide bond and a long aliphatic chain (ceramide), "
                  "and sugar rings; consistent with being a ganglioside.")

# Example usage:
if __name__ == "__main__":
    # Test using one of the provided SMILES strings.
    test_smiles = "O1[C@](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)C(O)[C@@H](O[C@@H]3CO)OC[C@H](NC(=O)CCCCCCCCCCCCCCC)[C@H](O)/C=C/CCCCCCCCCCCCC)(C[C@H](O)[C@@H](NC(=O)C)C1[C@H](O)[C@H](O)CO)C(O)=O"
    result, reason = is_ganglioside(test_smiles)
    print("Is ganglioside:", result)
    print("Reason:", reason)