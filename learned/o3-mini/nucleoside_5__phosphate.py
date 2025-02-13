"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate

Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which
            C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

Heuristic improvements over the prior version:
  • Identify a five-membered (furanose) ring with exactly one oxygen, four carbons.
  • In that ring, require exactly one exocyclic –CH2– group (CH2OH) that is linked via an oxygen
    to a “clean” phosphate chain (i.e. no extra carbon substituents on the phosphate(s)).
  • Verify that the sugar is also attached (through some atom of the ring) to an external
    aromatic heterocycle containing at least one nitrogen (and not itself a sugar ring).
  • Use a relatively relaxed molecular weight cutoff since valid nucleotides (mono-, di- or tri-phosphates)
    often exceed 600 Da but weed out very large molecules.
    
If the heuristic conditions are all met, the molecule is reported as a nucleoside 5'-phosphate.
Otherwise, the function returns False with a brief explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate from its SMILES string.
    
    A nucleoside 5'-phosphate is defined as a ribosyl or deoxyribosyl derivative of a
    purine or pyrimidine base in which the 5'-position of the sugar is phosphorylated.
    
    The heuristic tests applied are:
      1. Find a candidate sugar ring: a 5-membered ring with exactly one oxygen and four carbons.
      2. In that sugar, require exactly one exocyclic CH2 group (as seen in a CH2OH group).
         From that CH2, follow the connection (an oxygen) to verify a “clean” phosphate chain –
         that is, a phosphorus (or chain of phosphorus atoms) whose non-linking neighbors are only oxygens (or phosphorus).
      3. Check that the sugar ring connects (via any atom in the ring) to an external nucleobase:
         an external ring (ring size >= 5) that is aromatic, contains at least one nitrogen,
         and is not itself another sugar.
      4. As a safeguard, require that the overall molecular weight is not excessively high.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a nucleoside 5'-phosphate, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    
    # Use a moderate upper limit to weed out very large molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a typical nucleoside 5'-phosphate"
    
    ######## Helper functions ########
    def is_candidate_sugar_ring(ring):
        # Check ring length exactly 5 and count atoms
        if len(ring) != 5:
            return False
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        c_count   = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        return oxy_count == 1 and c_count == 4

    def is_CH2_group(atom):
        # Check if atom is carbon, sp3 hybridized, with exactly two attached hydrogens.
        if atom.GetAtomicNum() != 6:
            return False
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            return False
        # total num H (explicit+implicit)
        if atom.GetTotalNumHs() != 2:
            return False
        return True

    def verify_phosphate_chain(o_atom, visited_p=set()):
        """
        Given an oxygen atom o_atom that should be linking to a phosphate,
        check that at least one of its neighbors is a phosphorus that – aside from the linking oxygen –
        has only oxygen(s) or further phosphorus atoms as substituents.
        The function recurses on bridging phosphorus atoms.
        """
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 15:  # phosphorus
                if nbr.GetIdx() in visited_p:
                    continue
                visited_p.add(nbr.GetIdx())
                valid = True
                for p_neigh in nbr.GetNeighbors():
                    # Skip the oxygen that connected to it (o_atom)
                    if p_neigh.GetIdx() == o_atom.GetIdx():
                        continue
                    # Allow oxygen or phosphorus (for chained phosphates)
                    if p_neigh.GetAtomicNum() == 8:
                        continue
                    elif p_neigh.GetAtomicNum() == 15:
                        # recursively check the next phosphorus
                        if p_neigh.GetIdx() not in visited_p:
                            valid = valid and verify_phosphate_chain(p_neigh, visited_p)
                    else:
                        valid = False
                        break
                if valid:
                    return True
        return False

    def has_nucleobase(sugar_ring):
        # For each atom in the sugar ring, check if any neighbor (outside the ring)
        # is part of an aromatic ring (of size>=5) that contains at least one nitrogen.
        for idx in sugar_ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in sugar_ring:
                    continue
                # Check every ring that nbr belongs to.
                for ring in ring_info.AtomRings():
                    if nbr.GetIdx() in ring and len(ring) >= 5:
                        # Exclude rings that look like a sugar (i.e. 5-membered with 1 oxygen).
                        if is_candidate_sugar_ring(ring):
                            continue
                        # Check ring is aromatic and has at least one nitrogen.
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring) and \
                           any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in ring):
                            return True
        return False

    ######## End helper functions ########
    
    # Find candidate sugar rings.
    candidate_sugars = []
    for ring in ring_info.AtomRings():
        if is_candidate_sugar_ring(ring):
            candidate_sugars.append(ring)
    
    if not candidate_sugars:
        return False, "No 5-membered ring with 1 oxygen (expected for a ribofuranose) found"
    
    sugar_found = False
    reason_msg = ""
    
    # For each candidate sugar, look for the exocyclic CH2 that should lead (via O) to a phosphate chain.
    for sugar in candidate_sugars:
        phosphate_found = False
        exo_CH2_found = False
        
        # We expect exactly one exocyclic CH2 group (the 5'-CH2) attached to the ring.
        for idx in sugar:
            sugar_atom = mol.GetAtomWithIdx(idx)
            for nbr in sugar_atom.GetNeighbors():
                if nbr.GetIdx() in sugar:
                    continue
                # Look for a CH2 group. (Note: sometimes explicit hydrogens may be missing,
                # so we rely on GetTotalNumHs.)
                if is_CH2_group(nbr):
                    exo_CH2_found = True
                    # Now, from that CH2, check for an attached oxygen (should be –CH2OH)
                    found_O = False
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetIdx() == sugar_atom.GetIdx():
                            continue
                        if subnbr.GetAtomicNum() == 8:
                            # Check that this oxygen is linked to a (chain of) phosphate(s)
                            if verify_phosphate_chain(subnbr, set()):
                                found_O = True
                                phosphate_found = True
                                break
                    if found_O:
                        break
            if phosphate_found:
                break
        
        if not exo_CH2_found:
            reason_msg = "No exocyclic CH2 group (expected 5'-CH2 of the sugar) found"
            continue
        if not phosphate_found:
            reason_msg = "Exocyclic CH2 group did not lead to a proper phosphate chain"
            continue
        
        # Now check for attachment to a nucleobase.
        if not has_nucleobase(sugar):
            reason_msg = "No nucleobase (aromatic heterocycle with nitrogen) attached to the sugar was found"
            continue
        
        # If both phosphate and base are attached properly, we found our candidate.
        sugar_found = True
        break

    if sugar_found:
        return True, "Structure contains a nucleoside 5'-phosphate moiety with appropriate sugar, phosphate, and nucleobase"
    else:
        if reason_msg == "":
            reason_msg = "No candidate sugar ring scored positively for both phosphate and nucleobase attachments"
        return False, reason_msg

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O"  # adenosine 5'-monophosphate
    res, reason = is_nucleoside_5__phosphate(test_smiles)
    print(f"Test AMP: {res} // {reason}")
    
    # Some test examples from the list:
    ex_smiles = [
        "N1=C(C)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O",  # 6-azathymidine 5'-monophosphate (should be True)
        "CN1CN([C@@H]2O[C@H](COP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]2O)c2nc(N)[nH]c(=O)c12",  # dUTP (True)
        "C[C@@H]1O[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H](C[C@@H]2O)n2cc(C)c(=O)[nH]c2=O)",  # false positive example dTDP — not fitting our HEURISTICS
    ]
    
    for s in ex_smiles:
        r, msg = is_nucleoside_5__phosphate(s)
        print(f"SMILES: {s}\nResult: {r} // {msg}\n")