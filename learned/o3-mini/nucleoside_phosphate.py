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
  1. Identifying a candidate sugar ring (a 5-membered ring with 4 carbons and 1 oxygen).
  2. In addition, extending the sugar detection to include exocyclic hydroxyls (oxygen atoms
     attached to a ring carbon) that may be phosphorylated.
  3. Checking for a phosphorus atom (with neighbors including oxygen) attached to one of those
     exocyclic oxygens.
  4. Looking for a nucleobase by trying several SMARTS patterns (purine and pyrimidine cores).
Because of the molecular diversity, this remains a heuristic classification.
"""

from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    
    The algorithm expects to find:
       a) A candidate sugar ring (heuristically a furanose: 5 atoms, 4 carbons and 1 oxygen)
       b) A phosphate group (a phosphorus atom bonded via an oxygen to an exocyclic hydroxyl on the sugar)
       c) A nucleobase substructure (using several SMARTS for purine or pyrimidine rings).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a nucleoside phosphate, False otherwise.
        str: Explanation/reason for classification.
    """
    # Parse the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Identify candidate furanose rings (5-membered ring with 4 carbons and 1 oxygen).
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # List of sets of atom indices corresponding to candidate sugar rings.
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
    
    # 2. Check if a phosphate group is attached to the sugar moiety.
    # Instead of looking for phosphorus directly attached to atoms in the ring (which would fail for exocyclic attachments),
    # we look for an oxygen attached to a sugar ring carbon (exocyclic hydroxyl) that itself is bonded to a phosphorus atom.
    phosphate_attached = False
    # For each candidate sugar ring:
    for sugar in sugar_rings:
        # Look at each atom in the ring that is a carbon (possible attachment point for hydroxyls)
        for idx in sugar:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Check neighbors not in the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in sugar:
                    continue
                # We expect the hydroxyl oxygen to be exocyclic; it should be an oxygen.
                if nbr.GetAtomicNum() == 8:
                    # Now check if this oxygen (nbr) is connected to a phosphorus.
                    for o_nbr in nbr.GetNeighbors():
                        if o_nbr.GetIdx() == atom.GetIdx():
                            continue
                        if o_nbr.GetAtomicNum() == 15:
                            phosphate_attached = True
                            break
                    if phosphate_attached:
                        break
            if phosphate_attached:
                break
        if phosphate_attached:
            break
    if not phosphate_attached:
        return False, "No phosphate group found attached to a hydroxyl of the sugar moiety"
    
    # 3. Check for the presence of a nucleobase.
    # A list of SMARTS patterns covering typical purine and pyrimidine cores.
    nucleobase_smarts_list = [
        "n1cnc2ncnc12",           # common purine motif (adenine/guanine core)
        "c1nc2c(n1)cnc2",          # alternative purine pattern
        "n1cnc(=O)n1",            # cytosine-like pattern
        "O=c1[nH]c(=O)[nH]1",      # uracil-like pattern
        "C[C;H3]1=CN(C(=O)NC1=O)"   # pattern to pick up thymine-like structures (with methyl) (heuristic)
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
    
    return True, "Molecule contains a candidate sugar with a phosphate attached via an exocyclic hydroxyl and a nucleobase substructure."

# Example usage (you can run these tests when invoking this script):
if __name__ == "__main__":
    test_molecules = {
        "GTP" : "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O",
        "UDP" : "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O",
        "N(4)-acetylcytidine 5'-monophosphate": "CC(=O)Nc1ccn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c(=O)n1"
    }
    
    for name, smiles in test_molecules.items():
        result, reason = is_nucleoside_phosphate(smiles)
        print(f"{name}: {result}\n  Reason: {reason}\n")