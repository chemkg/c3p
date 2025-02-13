"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate

Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which
            C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.
            
Heuristic Approach:
  1. Parse the SMILES and get ring info.
  2. Search for a candidate “sugar” ring. In a nucleoside sugar (ribose or deoxyribose),
     the cyclic (furanose) ring is expected to have 5 atoms, exactly one of which
     is oxygen and the other four carbons.
  3. Verify that at least one atom of that sugar is attached outside the ring to:
       a. a substituent leading to a phosphate group (or a phosphate chain) that is “pure”
          (i.e. all non-connection atoms at the phosphate(s) are oxygen – no extra carbon appendage),
       b. and a nucleobase (i.e. an external aromatic ring that contains at least one nitrogen).
  4. Finally, as a safeguard to weed out large molecules (e.g. CoA derivatives) that merely contain
     a nucleoside 5'-phosphate substructure, we check that the molecular weight is not too high.
     
If all criteria are met then the structure is classified as a nucleoside 5'-phosphate.
If there is an issue with the substructure or the molecule is too large, the function returns False.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate from its SMILES string.
    
    A nucleoside 5'-phosphate is defined as a ribosyl or deoxyribosyl derivative of a
    purine or pyrimidine base in which the 5'-position of the sugar is phosphorylated
    (mono-, di-, tri- or tetra-phosphorylated).
    
    We apply the following heuristic tests:
      1. Identify a 5-membered ring with exactly one oxygen and four carbons (a furanose ring).
      2. In that candidate sugar, check for an exocyclic substituent (expected to be the 5'-CH2)
         that is attached via an oxygen to a phosphate group (or chain of phosphate groups).
           • For the attached phosphate, require that every neighbor (other than the linking oxygen)
             is an oxygen or another phosphate (indicating an “unadorned” phosphate chain).
      3. Check that the sugar ring is connected to a nucleobase: at least one sugar atom must have 
         a neighbor (outside the ring) that belongs to an external ring (size>=5) that is aromatic
         and contains at least one nitrogen.
      4. Finally, if the overall molecular weight is unusually high (e.g. >600 Da) we assume this is
         part of a larger molecule rather than a pure nucleoside 5'-phosphate.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a nucleoside 5'-phosphate, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    
    # Check overall molecular weight: typical nucleoside mono/diphosphates are below ~600 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 600:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a typical nucleoside 5'-phosphate"
    
    sugar_candidate_found = False
    phosphate_attached = False
    nucleobase_attached = False

    # Step 1: find candidate sugar ring (a 5-membered ring with exactly one oxygen and four carbons)
    sugar_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        c_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if oxy_count == 1 and c_count == 4:
            sugar_rings.append(ring)

    if not sugar_rings:
        return False, "No 5-membered ring with 1 oxygen (expected for a ribofuranose) found"
    
    # For each candidate sugar ring, try to verify the required external connections.
    for ring in sugar_rings:
        # Step 2a: Look for a phosphate connection.
        # We expect that one of the sugar atoms is connected (exocyclically) to a CH2 group,
        # which in turn is linked via an oxygen to a phosphate.
        this_ring_has_phosphate = False
        for idx in ring:
            sugar_atom = mol.GetAtomWithIdx(idx)
            for nbr in sugar_atom.GetNeighbors():
                # Consider neighbors outside the ring
                if nbr.GetIdx() in ring:
                    continue
                # Expect a carbon (typically primary CH2) emanating from the sugar
                if nbr.GetAtomicNum() != 6 or nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                # Now check if this exocyclic carbon has a neighbor that is an oxygen attached to a phosphate
                for subnbr in nbr.GetNeighbors():
                    # We want an oxygen that is not the sugar carbon (nbr)
                    if subnbr.GetAtomicNum() != 8:
                        continue
                    # Search neighbors of that oxygen for a phosphorus atom
                    for p_neigh in subnbr.GetNeighbors():
                        if p_neigh.GetAtomicNum() == 15:
                            # Check that—aside from the linking oxygen—every other neighbor of the phosphorus (or chain)
                            # is oxygen or another phosphorus (i.e. no additional carbon, which would indicate an extra
                            # acyl or conjugated group).
                            valid_phosphate = True
                            for pn in p_neigh.GetNeighbors():
                                if pn.GetIdx() == subnbr.GetIdx():
                                    continue
                                if pn.GetAtomicNum() not in (8, 15):
                                    valid_phosphate = False
                                    break
                            if valid_phosphate:
                                this_ring_has_phosphate = True
                                break
                    if this_ring_has_phosphate:
                        break
                if this_ring_has_phosphate:
                    break
            if this_ring_has_phosphate:
                break
        
        # Step 2b: Look for attachment to a nucleobase.
        # We require that at least one atom of the sugar ring is directly connected to an external ring
        # (of size >= 5) that is aromatic and contains at least one nitrogen.
        this_ring_has_base = False
        for idx in ring:
            atom_in_ring = mol.GetAtomWithIdx(idx)
            for nbr in atom_in_ring.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # For each neighbor outside the sugar, check if it appears in any sufficiently large aromatic ring.
                for other_ring in ring_info.AtomRings():
                    if nbr.GetIdx() in other_ring and len(other_ring) >= 5:
                        # Check all atoms in this ring are aromatic and at least one is a nitrogen.
                        if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in other_ring) and \
                           any(mol.GetAtomWithIdx(a).GetAtomicNum() == 7 for a in other_ring):
                            this_ring_has_base = True
                            break
                if this_ring_has_base:
                    break
            if this_ring_has_base:
                break
        
        if this_ring_has_phosphate and this_ring_has_base:
            sugar_candidate_found = True
            phosphate_attached = True
            nucleobase_attached = True
            break  # our candidate sugar meets everything

    if not sugar_candidate_found:
        msg = "Did not find a candidate sugar ring with both a terminal phosphate and attached nucleobase"
        if not phosphate_attached:
            msg = "No proper phosphate (or phosphate chain) attached to a sugar exocyclic group was found"
        elif not nucleobase_attached:
            msg = "No nucleobase (aromatic heterocycle with nitrogen) attached to the sugar was found"
        return False, msg

    return True, "Structure contains a nucleoside 5'-phosphate moiety with appropriate sugar, phosphate, and nucleobase"

# Example usage (for testing):
if __name__ == "__main__":
    # A known nucleoside 5'-monophosphate (e.g., adenosine 5'-monophosphate)
    test_smiles = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O"
    res, reason = is_nucleoside_5__phosphate(test_smiles)
    print(f"Test AMP: {res} // {reason}")
    
    # A false positive candidate (a CoA derivative fragment - note: simplified SMILES)
    false_smiles = "C[C@H](CCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](O)[C@@H]1OP(O)(O)=O)"
    res, reason = is_nucleoside_5__phosphate(false_smiles)
    print(f"Test false candidate: {res} // {reason}")