"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
#!/usr/bin/env python3
"""
Classifies: N-acyl-L-alpha-amino acid
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An N-acyl-L-alpha-amino acid should have one “free” L-alpha-amino acid core 
(i.e. a chiral alpha‑carbon with an attached free carboxyl group, not in a peptide bond)
and at least one nitrogen (which may be the alpha amine or a side-chain amine) that is acylated.
Acylation is defined by a neighboring carbon (other than the carboxyl carbon) that bears a C=O bond
and is not part of a peptide bond.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    The approach is:
      1. Find all potential L-alpha amino acid cores, that is, chiral carbons (with '@')
         with an attached carboxyl group (C(=O)O or C(=O)[O-]). (We allow both neutral and deprotonated groups.)
         Then, try to filter out cores that are clearly part of a peptide chain by searching for
         an amide connection to another L-alpha core.
      2. From the remaining core(s), choose one candidate (if more than one remain, pick the one that
         has an acylated nitrogen by further analysis).
      3. Check for an acylated nitrogen (within a topological distance up to 5 bonds from the core)
         that is not part of a peptide linkage. A nitrogen is acylated if it is bonded to a carbon
         (other than the carboxyl carbon on the candidate core) that carries a double-bonded oxygen.
         
    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      bool: True if the molecule is classified as an N-acyl-L-alpha-amino acid, False otherwise.
      str: Explanation of the reasoning.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Identify potential L-alpha cores.
    # We use several SMARTS to allow both neutral and deprotonated carboxyl groups.
    core_smarts = [
        "[C@H](C(=O)O)",
        "[C@H](C(=O)[O-])",
        "[C@@H](C(=O)O)",
        "[C@@H](C(=O)[O-])"
    ]
    candidate_cores = []  # Each entry: dict with keys: alpha_idx and acid_idx.
    for smarts in core_smarts:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            # In our SMARTS the first atom is the chiral alpha carbon.
            alpha_idx = match[0]
            # Identify the carboxyl carbon among neighbors of the alpha atom.
            alpha_atom = mol.GetAtomWithIdx(alpha_idx)
            acid_idx = None
            for nbr in alpha_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:  # carbon
                    # Look for a double-bonded oxygen on this neighbor.
                    for nn in nbr.GetNeighbors():
                        if nn.GetAtomicNum() == 8:
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                            if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                                acid_idx = nbr.GetIdx()
                                break
                    if acid_idx is not None:
                        break
            if acid_idx is not None:
                candidate_cores.append({'alpha_idx': alpha_idx, 'acid_idx': acid_idx})
    
    if not candidate_cores:
        return False, "No L-alpha-amino acid core (chiral carbon with attached carboxyl) detected"
    
    # Step 1b. Filter out cores that appear to be part of a peptide chain.
    # We define a peptide bond if the carboxyl carbon (of one core) is connected
    # to a nitrogen which in turn is attached to another potential core’s alpha carbon.
    def is_core_peptide_link(core, other_core):
        # Check if the acid carbon of core is bonded to a nitrogen that is also bonded to other_core's alpha.
        acid_atom = mol.GetAtomWithIdx(core['acid_idx'])
        for nbr in acid_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 7:  # nitrogen
                # Check if this nitrogen is bonded to the other core's alpha carbon.
                for nn in nbr.GetNeighbors():
                    if nn.GetIdx() == other_core['alpha_idx']:
                        return True
        return False

    free_cores = []
    for i, core in enumerate(candidate_cores):
        peptide_linked = False
        for j, other in enumerate(candidate_cores):
            if i == j:
                continue
            if is_core_peptide_link(core, other):
                peptide_linked = True
                break
        # Allow cores that are not peptide-linked or that are connected only via a side chain (see below).
        if not peptide_linked:
            free_cores.append(core)
    
    # If no free core is found, we assume the molecule is a peptide or not a simple amino acid.
    if not free_cores:
        return False, "More than one L-alpha-amino acid core found in a peptide-like arrangement"
    
    # If more than one free core remains, try to choose one which (upon later checks) shows acylation.
    # (We loop over free cores in the next step.)
    
    # Helper: Check if a nitrogen atom is acylated.
    def is_acylated_nitrogen(n_atom, exclude_acid_idx):
        # A nitrogen is considered acylated if any of its neighboring carbons (other than the exclude acid)
        # has a double-bonded oxygen and is not clearly peptide-bonded.
        for nbr in n_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != exclude_acid_idx:
                has_carbonyl = False
                for nn in nbr.GetNeighbors():
                    if nn.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                            has_carbonyl = True
                            break
                if has_carbonyl:
                    # Check: if this acyl carbon is bonded to any nitrogen that is connected to a free core
                    # (other than via the current n_atom) then it may be a peptide bond.
                    peptide_flag = False
                    for nn in nbr.GetNeighbors():
                        if nn.GetAtomicNum() == 7 and nn.GetIdx() != n_atom.GetIdx():
                            for fc in free_cores:
                                if fc['alpha_idx'] == nn.GetIdx():
                                    peptide_flag = True
                                    break
                        if peptide_flag:
                            break
                    if not peptide_flag:
                        return True
        return False

    # Now, for each free core candidate, look for a nitrogen (in the vicinity
    # of the core – we search all atoms within a maximum path length) that is acylated.
    max_distance = 5
    from collections import deque

    def find_acylated_nitrogen_from_core(core):
        start_idx = core['alpha_idx']
        visited = {start_idx}
        queue = deque([(start_idx, 0)])
        while queue:
            current, dist = queue.popleft()
            if 0 < dist <= max_distance:
                atom = mol.GetAtomWithIdx(current)
                if atom.GetAtomicNum() == 7:
                    # If this nitrogen is directly attached to the alpha carbon, then
                    # exclude the acid carbon from its acylation check.
                    if is_acylated_nitrogen(atom, core['acid_idx']):
                        return (atom.GetIdx(), dist)
            if dist < max_distance:
                atom = mol.GetAtomWithIdx(current)
                for nbr in atom.GetNeighbors():
                    idx = nbr.GetIdx()
                    if idx not in visited:
                        visited.add(idx)
                        queue.append((idx, dist+1))
        return None

    # Loop over free cores. If a core gives us an acylated nitrogen, we accept and return True.
    for core in free_cores:
        res = find_acylated_nitrogen_from_core(core)
        if res is not None:
            n_idx, dist = res
            return True, (
                "Molecule has an L-alpha-amino acid core (alpha carbon index %d with attached carboxyl carbon %d) "
                "and an acylated nitrogen (atom index %d within %d bonds) that is not part of a peptide bond."
                % (core['alpha_idx'], core['acid_idx'], n_idx, max_distance)
            )

    return False, "No qualifying N-acyl substituent found attached to an L-alpha-amino acid core."

# Example testing:
if __name__ == '__main__':
    # You can test with one of the provided SMILES strings.
    test_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"  # N-acetyl-L-aspartic acid
    result, reason = is_N_acyl_L_alpha_amino_acid(test_smiles)
    print(result, reason)