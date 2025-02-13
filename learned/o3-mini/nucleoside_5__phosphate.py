"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate

Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which
            C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

Heuristic improvements over the prior version:
  • Identify a candidate sugar ring: a 5-membered ring with exactly one oxygen and four carbons.
  • Within that ring, look for one or more exocyclic CH2 groups (expected for the 5'-CH2OH).
    For at least one such CH2, require that one of its oxygen substituents is directly linked
    to a “clean” phosphate chain (i.e. the connected phosphorus – and any chained phosphorus –
    must have only oxygen or phosphorus as neighbors besides the linking oxygen).
  • Check that the sugar ring is attached (via any atom of the ring) to an external nucleobase.
    Here the nucleobase is defined as an aromatic ring (size >= 5) containing at least one nitrogen;
    we relax the requirement slightly by allowing a neighbor that is aromatic and (ideally) a nitrogen.
  • Use a moderate molecular weight cutoff to weed out excessively large molecules.

If all conditions are met, the molecule is classified as nucleoside 5'-phosphate; otherwise,
the function returns False with a brief explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate from its SMILES string.
    
    A nucleoside 5'-phosphate is defined as a ribosyl or deoxyribosyl derivative of a
    purine or pyrimidine base where the 5'-position of the sugar is phosphorylated.

    The heuristic tests are:
      1. Look for a candidate sugar: a 5-membered ring containing exactly one oxygen and four carbons.
      2. In the sugar, find at least one exocyclic CH2 group (expected as the 5'-CH2).
         For at least one such CH2 group, check that one oxygen substituent is directly linked
         to a phosphate chain. In our phosphate chain check, we require that the phosphorus
         (or connected phosphorus atoms) has only oxygen or phosphorus as its other substituents.
      3. Verify that at least one atom of the sugar ring is bound to a nucleobase. In our test,
         we require that a neighbor (outside the sugar) belongs to an aromatic ring (of size ≥5)
         that contains at least one nitrogen and is not itself a sugar ring.
      4. Reject overly heavy molecules.
      
    Args:
      smiles (str): SMILES string
    
    Returns:
      bool: True if classified as nucleoside 5'-phosphate, else False.
      str: Explanation message.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    
    # Reject very large molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a typical nucleoside 5'-phosphate"
    
    ############ Helper functions ############
    
    def is_candidate_sugar_ring(ring):
        # A candidate sugar ring must have exactly 5 atoms
        if len(ring) != 5:
            return False
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        c_count   = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        return oxy_count == 1 and c_count == 4

    def is_CH2_group(atom):
        # Returns True if the atom is carbon (atomic # 6), sp3 hybridized and has exactly 2 hydrogens.
        if atom.GetAtomicNum() != 6:
            return False
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            return False
        if atom.GetTotalNumHs() != 2:
            return False
        return True
    
    def verify_phosphate_chain(o_atom, visited_p=set()):
        """
        Given an oxygen atom o_atom (coming from a candidate CH2 group),
        check if one of its neighbors is a phosphorus that is "clean":
        all other substituents on that phosphorus (or connected phosphorus atoms
        recursively) must be either oxygen (atomic number 8) or phosphorus (15).
        """
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 15:  # phosphorus
                if nbr.GetIdx() in visited_p:
                    continue
                visited_p.add(nbr.GetIdx())
                valid = True
                for p_neigh in nbr.GetNeighbors():
                    # Skip the oxygen that connected to it
                    if p_neigh.GetIdx() == o_atom.GetIdx():
                        continue
                    # Only allow oxygen or phosphorus (chained phosphates)
                    if p_neigh.GetAtomicNum() not in (8, 15):
                        valid = False
                        break
                    # If neighbor is phosphorus, check recursively
                    if p_neigh.GetAtomicNum() == 15 and p_neigh.GetIdx() not in visited_p:
                        valid = valid and verify_phosphate_chain(p_neigh, visited_p)
                if valid:
                    return True
        return False
    
    def has_nucleobase(sugar_ring):
        """
        For each atom in the candidate sugar ring, examine its neighbors (outside the ring)
        to see if it is connected to an aromatic (ring size ≥ 5) system containing at least one nitrogen.
        We relax the condition by accepting a neighbor that is aromatic and, when examined in the context
        of its ring, has at least one nitrogen.
        """
        for idx in sugar_ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in sugar_ring:
                    continue
                # If neighbor is aromatic and is a nitrogen, return True immediately.
                if nbr.GetAtomicNum() == 7 and nbr.GetIsAromatic():
                    return True
                # Otherwise look in every ring that this neighbor is part of.
                for ring in ring_info.AtomRings():
                    if nbr.GetIdx() in ring and len(ring) >= 5:
                        # Skip if the ring is likely another sugar ring.
                        if is_candidate_sugar_ring(ring):
                            continue
                        # Check that at least one atom in this ring is a nitrogen
                        if any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in ring):
                            return True
        return False
    
    ############ End helper functions ############
    
    # Find candidate sugar rings.
    candidate_sugars = []
    for ring in ring_info.AtomRings():
        if is_candidate_sugar_ring(ring):
            candidate_sugars.append(ring)
    
    if not candidate_sugars:
        return False, "No 5-membered ring with 1 oxygen (expected for a ribofuranose) found"
    
    sugar_valid = False
    msg_reason = ""
    
    # Look for at least one candidate sugar that has both a properly linked phosphate chain and an attached nucleobase.
    for sugar in candidate_sugars:
        phosphate_found = False
        exo_CH2_found = False
        
        # Look among the neighbors of sugar atoms for an exocyclic CH2 group.
        for idx in sugar:
            sugar_atom = mol.GetAtomWithIdx(idx)
            for nbr in sugar_atom.GetNeighbors():
                if nbr.GetIdx() in sugar:
                    continue  # already in sugar ring
                if is_CH2_group(nbr):
                    exo_CH2_found = True
                    # Check that one of its oxygen neighbors leads cleanly to a phosphorus.
                    valid_link = False
                    for subnbr in nbr.GetNeighbors():
                        # Exclude the atom through which CH2 attached to the sugar.
                        if subnbr.GetIdx() == sugar_atom.GetIdx():
                            continue
                        if subnbr.GetAtomicNum() == 8:
                            # Ensure that this oxygen leads directly to a phosphate group.
                            if verify_phosphate_chain(subnbr, set()):
                                valid_link = True
                                break
                    if valid_link:
                        phosphate_found = True
                        break
            if phosphate_found:
                break
        
        if not exo_CH2_found:
            msg_reason = "No exocyclic CH2 group (expected 5'-CH2 of the sugar) found"
            continue
        if not phosphate_found:
            msg_reason = "Exocyclic CH2 group did not lead directly to a proper phosphate group"
            continue
        
        if not has_nucleobase(sugar):
            msg_reason = "No nucleobase (aromatic heterocycle with nitrogen) attached to the sugar was found"
            continue
        
        # All tests passed for this sugar.
        sugar_valid = True
        break

    if sugar_valid:
        return True, "Structure contains a nucleoside 5'-phosphate moiety with appropriate sugar, phosphate, and nucleobase"
    else:
        if msg_reason == "":
            msg_reason = "No candidate sugar ring met the requirements for nucleotide structure"
        return False, msg_reason

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O"  # adenosine 5'-monophosphate
    res, reason = is_nucleoside_5__phosphate(test_smiles)
    print(f"Test AMP: {res} // {reason}")
    
    # Some test examples from the list:
    examples = [
        ("6-azathymidine 5'-monophosphate", "N1=C(C)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O"),
        ("dUTP", "O[C@H]1C[C@@H](O[C@@H]1COP(O)(=O)OP(O)(=O)OP(O)(O)=O)N1C=CC(=O)NC1=O"),
        ("5-carboxy-2'-deoxyuridine 5'-monophosphate", "[C@@H]1(N2C(NC(=O)C(=C2)C(=O)O)=O)O[C@H](COP(O)(=O)O)[C@H](C1)O"),
    ]
    for name, s in examples:
        res, msg = is_nucleoside_5__phosphate(s)
        print(f"SMILES for {name}: {s}\nResult: {res} // {msg}\n")