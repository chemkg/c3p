"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
#!/usr/bin/env python
"""
Classifies: nucleoside phosphate

A nucleoside phosphate is defined as: 
  "A nucleobase-containing molecular entity that is a nucleoside in which one or more of the 
   sugar hydroxy groups has been converted into a mono- or poly-phosphate. The term includes 
   both nucleotides and non-nucleotide nucleoside phosphates."

This implementation improves upon earlier attempts by:
  1. Identifying a candidate sugar ring (a furanose: 5-membered ring with 4 carbons and 1 oxygen).
  2. Expanding the sugar region to include atoms directly attached to the ring (such as the exocyclic 5'-CH2).
  3. Looking for a phosphorus atom whose neighbor oxygen is directly attached to the sugar region.
  4. Searching for a nucleobase substructure using several SMARTS patterns for purine and pyrimidine cores.
Because of the diversity of nucleoside phosphates, this remains a heuristic classification.
"""

from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.

    The algorithm uses the following steps:
       a) Identify a candidate sugar ring: a 5-membered ring with exactly one oxygen and four carbons.
       b) Expand the sugar region by including atoms directly bonded to the sugar ring.
       c) Check that at least one phosphorus (P) is attached via an oxygen to an atom in the sugar region.
       d) Verify that a nucleobase-like substructure is present (using multiple SMARTS for purine/pyrimidine cores).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple with the boolean classification and an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Identify candidate sugar rings.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # each is a set of atom indices that belong to the candidate sugar ring.
    for ring in ring_info:
        if len(ring) == 5:
            o_count = 0
            c_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    o_count += 1
                elif atom.GetAtomicNum() == 6:
                    c_count += 1
            if o_count == 1 and c_count == 4:
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No candidate sugar ring (5-membered ring with 4 carbons and 1 oxygen) detected"

    # Step 2: Expand the sugar region.
    # For each candidate sugar ring found, add the atoms directly bonded to atoms in the ring.
    sugar_region = set()
    for ring in sugar_rings:
        sugar_region |= ring  # include all atoms in the ring
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                sugar_region.add(nbr.GetIdx())
    # (This expansion should capture exocyclic groups such as a 5'-CH2)

    # Step 3: Check for phosphate group attached to the sugar region.
    phosphate_found = False
    # Look for any phosphorus (atomic num 15)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            # For each phosphorus, look at oxygen neighbors.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen is connected to an atom in the sugar region.
                    for o_nbr in nbr.GetNeighbors():
                        if o_nbr.GetIdx() in sugar_region and o_nbr.GetAtomicNum() in (6, 8):
                            phosphate_found = True
                            break
                    if phosphate_found:
                        break
            if phosphate_found:
                break
    if not phosphate_found:
        return False, "No phosphate group found attached (via an oxygen) to the sugar region"

    # Step 4: Check for the presence of a nucleobase substructure.
    # We include several SMARTS for common purine and pyrimidine patterns.
    nucleobase_smarts_list = [
        "n1cnc2ncnc12",              # common purine (adenine/guanine) motif
        "c1ncnc2ncnc12",             # alternative purine motif
        "Nc1nc2n(cnc2c(=O)[nH]1)",    # guanine-like motif
        "O=c1[nH]cnc1=O",            # pyrimidine (uracil-like)
        "O=c1nc[nH]c1"               # cytosine-like
    ]
    nucleobase_found = False
    for smarts in nucleobase_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No nucleobase substructure found (expected purine or pyrimidine ring)"
    
    return True, "Molecule contains a candidate sugar with a phosphate attached via an oxygen, and a nucleobase substructure."

# Optional testing code:
if __name__ == "__main__":
    test_molecules = {
        "GTP": "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O",
        "UDP": "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O",
        "N(4)-acetylcytidine 5'-monophosphate": "CC(=O)Nc1ccn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c(=O)n1"
    }
    for name, smi in test_molecules.items():
        result, reason = is_nucleoside_phosphate(smi)
        print(f"{name}: {result}")
        print(f"  Reason: {reason}")