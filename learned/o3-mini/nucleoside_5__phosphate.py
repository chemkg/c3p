"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate
Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which
            C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.
            
This code uses heuristic SMARTS and ring analyses:
  1. It identifies a phosphate (P with atomic num 15) attached through an oxygen (atomic num 8)
     to a carbon.
  2. It then checks that this carbon is a member of a 5-membered ring that contains exactly one oxygen
     (consistent with a ribofuranose or deoxyribofuranose ring).
  3. It then confirms that at least one atom in that sugar ring is attached to an external aromatic ring
     that contains nitrogen atoms (indicative of a purine/pyrimidine nucleobase).
If these requirements are met, the molecule is classified as a nucleoside 5'-phosphate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate from its SMILES string.
    
    A nucleoside 5'-phosphate is defined as a ribosyl or deoxyribosyl derivative of a
    pyrimidine or purine base in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.
    
    Steps performed:
      1. Parse the SMILES string.
      2. Look for a phosphate (P) atom attached via an oxygen to a candidate sugar carbon.
      3. Among rings containing that sugar carbon, look for a 5-membered ring with exactly one oxygen (sugar ring).
      4. Check that the sugar ring is connected to an external aromatic heterocycle (the nucleobase).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a nucleoside 5'-phosphate, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    
    # Flag to mark if the required fingerprint was found
    nucleoside_found = False
    
    # Loop over all atoms to find a phosphorus (P) atom (atomic num 15)
    for p in mol.GetAtoms():
        if p.GetAtomicNum() != 15:
            continue  # skip if not phosphorus
        
        # Look at neighbors for oxygen atoms (atomic num 8)
        for o in p.GetNeighbors():
            if o.GetAtomicNum() != 8:
                continue
            # Now, from the oxygen (the bridging atom), find a carbon neighbor (candidate sugar carbon)
            for c in o.GetNeighbors():
                if c.GetAtomicNum() != 6:
                    continue
                # Ensure that this carbon is not the phosphorus itself
                if c.GetIdx() == p.GetIdx():
                    continue
                # Optionally check for sp3 hybridization (sugar carbons are sp3)
                if c.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                
                # Now check if this carbon (candidate for C-5) is part of a 5-membered ring 
                # that contains exactly one oxygen (expected for ribose/deoxyribose)
                candidate_ring = None
                for ring in ring_info.AtomRings():
                    if c.GetIdx() in ring and len(ring) == 5:
                        # Count oxygen atoms in the ring
                        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
                        if oxy_count == 1:
                            candidate_ring = ring
                            break
                if candidate_ring is None:
                    # No suitable sugar ring found that contains this candidate carbon
                    continue
                
                # Next, check that this sugar ring is connected to a nucleobase.
                # We do so by checking if any atom in the sugar ring has a neighbor not in the ring that is part 
                # of an external ring of size >=5, is aromatic, and contains at least one nitrogen.
                base_found = False
                for idx in candidate_ring:
                    sugar_atom = mol.GetAtomWithIdx(idx)
                    for nbr in sugar_atom.GetNeighbors():
                        if nbr.GetIdx() in candidate_ring:
                            continue  # neighbor is within the sugar ring itself
                        # Check if neighbor belongs to an external ring (candidate nucleobase ring)
                        for ring2 in ring_info.AtomRings():
                            if nbr.GetIdx() in ring2 and len(ring2) >= 5:
                                # Check if all atoms in ring2 are aromatic and that at least one is nitrogen
                                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring2):
                                    if any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in ring2):
                                        base_found = True
                                        break
                        if base_found:
                            break
                    if base_found:
                        break
                if not base_found:
                    continue  # sugar ring found but no connected nucleobase
                
                # If we reach this point then we have:
                # - A phosphate group (P) attached via an oxygen to a sugar carbon.
                # - A sugar ring (5-membered with one O) that is connected to an aromatic heterocycle.
                nucleoside_found = True
                break  # no need to search further for this P atom
            if nucleoside_found:
                break
        if nucleoside_found:
            break

    if not nucleoside_found:
        return False, "Molecule does not contain the expected nucleoside 5'-phosphate substructure"
    
    # Optional additional checks (molecular weight, rotatable bonds, etc.) can be added here if desired.
    
    return True, "Structure contains a nucleoside 5'-phosphate moiety (sugar-phosphate linkage and nucleobase detected)"