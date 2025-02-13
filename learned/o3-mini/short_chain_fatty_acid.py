"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: short-chain fatty acid

Definition:
  A short-chain fatty acid (SCFA) is an acyclic aliphatic monocarboxylic acid 
  that contains exactly one carboxyl group (matched as [CX3](=O)[OX1H]) and whose 
  longest continuous carbon chain (starting at the acid carbon) is between 3 and 5 carbons. 
  In addition, aside from the acid group, any substituent attached to a carbon in the parent chain 
  must be “small” – here we allow only a one‐atom branch (a methyl) or a branch that is a –CH2OH 
  (or an –OH directly attached).  (Thus formic and acetic acids are not considered SCFAs.)
  
Note: This is one interpretation to “improve” on the previous attempt.
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid according to the following criteria:
      1. Molecule must be acyclic.
      2. It must contain exactly one carboxyl group ([CX3](=O)[OX1H]).
      3. The longest chain of carbons (starting from the acid carbon) must have 3, 4, or 5 carbons.
      4. Any branch (i.e. substituent off the parent chain, except for the carboxyl oxygens)
         must be very small: allowed substituents are either:
           - a methyl group (a single carbon with no non-H neighbors beyond the attachment),
           - a –CH2OH group (i.e. one carbon that carries one –OH) or
           - an –OH group directly.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if molecule is a short-chain fatty acid per the above, False otherwise.
        str: Explanation for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Must be acyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and is not purely aliphatic"
    
    # Criterion 2: Must have exactly one carboxyl group.
    acid_query = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    acid_matches = mol.GetSubstructMatches(acid_query)
    if len(acid_matches) != 1:
        return False, "Molecule does not have exactly one carboxyl group"
    
    acid_match = acid_matches[0]  # (acid carbon, acid oxygen)
    acid_carbon = acid_match[0]
    acid_oxygens = set(acid_match[1:])  # carboxyl oxygens to be ignored on acid carbon

    # Build a graph over carbon atoms (only consider C–C bonds)
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_graph = {idx: [] for idx in carbon_atoms}
    for idx in carbon_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                nidx = nbr.GetIdx()
                carbon_graph[idx].append(nidx)

    if acid_carbon not in carbon_graph:
        return False, "Acid carbon not found in carbon graph"

    # We now enumerate all simple paths (no cycles) starting at acid carbon.
    # (Because the molecule is acyclic, the carbon graph is a tree.)
    candidate_paths = []
    def dfs_path(current, path):
        # record current path
        candidate_paths.append(path.copy())
        for nbr in carbon_graph[current]:
            if nbr not in path:
                path.append(nbr)
                dfs_path(nbr, path)
                path.pop()
    dfs_path(acid_carbon, [acid_carbon])
    
    if not candidate_paths:
        return False, "No carbon chain found from acid carbon"

    # We now select candidate parent chains whose length (number of carbons) is between 3 and 5.
    valid_parent_chains = []
    max_length = 0
    for path in candidate_paths:
        # length is the number of carbons in this simple path.
        chain_len = len(path)
        if chain_len > max_length:
            max_length = chain_len
        if 3 <= chain_len <= 5:
            valid_parent_chains.append(path)
    if not valid_parent_chains:
        return False, f"Longest carbon chain has {max_length} carbons, but valid SCFAs require 3-5 carbons in the chain"
    
    # A helper: check branches. For each carbon in the candidate chain,
    # check all neighbors that are not in the chain (and, for the acid carbon, ignore the acid oxygens)
    # Allowed branch types: 
    #   - If atom is carbon: it must be “methyl” (i.e. when looking at that branch, 
    #         it has no additional heavy-atom neighbors beyond the attachment) 
    #         OR it may be a CH2 that bears one oxygen substituent (–CH2OH).
    #   - If atom is oxygen: allowed only if it is -OH (i.e. attached only to H aside from the connection).
    def branch_allowed(nbr):
        # nbr is an atom (the neighbor off a backbone carbon)
        atomic_num = nbr.GetAtomicNum()
        if atomic_num == 6:
            # Examine neighbors of this branch atom other than the backbone connection.
            heavy_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() not in (1)]
            # We expect only one heavy neighbor if it is a simple methyl.
            if len(heavy_neighbors) == 1:
                return True  # a methyl branch
            elif len(heavy_neighbors) == 2:
                # Possibly a –CH2OH group.
                # The branch atom (should be CH2) along with exactly one additional heavy neighbor.
                others = [n for n in heavy_neighbors if n.GetIdx() not in branch_context]
                if len(others) == 1 and others[0].GetAtomicNum() == 8:
                    oxy = others[0]
                    # Ensure oxygen is -OH (its only heavy neighbor is this branch atom)
                    oxy_heavy = [t for t in oxy.GetNeighbors() if t.GetAtomicNum() not in (1) and t.GetIdx() != nbr.GetIdx()]
                    if len(oxy_heavy) == 0:
                        return True
            return False  # larger substituent not allowed
        elif atomic_num == 8:
            # If an oxygen is directly attached to the backbone: allow only if it is -OH.
            # Check that (ignoring H) oxygen has no other heavy neighbor.
            heavy_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() not in (1)]
            if len(heavy_neighbors) == 1:
                return True
            return False
        else:
            # Other atoms are not allowed in substituents.
            return False

    # Now check each candidate parent chain.
    # We set a flag if any candidate chain passes the branch tests.
    candidate_ok = False
    candidate_reason = ""
    for path in valid_parent_chains:
        backbone_set = set(path)
        ok = True
        # For each carbon in this chain, check its neighbors that are not in the chain.
        for idx in path:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nidx = nbr.GetIdx()
                # At the acid carbon, ignore carboxyl oxygens.
                if idx == acid_carbon and nidx in acid_oxygens:
                    continue
                # If neighbor is in the backbone, skip.
                if nidx in backbone_set:
                    continue
                # For each branch neighbor, check allowed substituent rules.
                # To help with look‐up in branch_allowed if more than one branch is connected,
                # we use a context (here we only want to check the immediate neighbor and not confuse it if it belongs to more than one branch).
                branch_context = backbone_set  # so that we only consider that neighbor outside the backbone.
                if not branch_allowed(nbr):
                    ok = False
                    candidate_reason = f"Branch on backbone carbon {idx} is not an allowed substituent"
                    break
            if not ok:
                break
        if ok:
            candidate_ok = True
            chosen_len = len(path)
            break

    if candidate_ok:
        return True, f"Molecule is an acyclic monocarboxylic acid with a {chosen_len}-carbon parent chain and only allowed substituents"
    else:
        return False, candidate_reason or "No valid parent chain with allowed substituents found"


# (Optional) Testing block
if __name__ == "__main__":
    test_examples = [
        # True positives
        ("CCCCC(O)=O", "valeric acid"),
        ("CC(C)=CC(O)=O", "3-methylbut-2-enoic acid"),
        ("CC(C)C(C)C(O)=O", "2,3-dimethylbutyric acid"),
        ("OC(=O)CC=C(C)C", "4-methylpent-3-enoic acid"),
        ("CC[C@@H](C)C(O)=O", "(R)-2-methylbutyric acid"),
        ("CCCC(C)C(O)=O", "2-methylvaleric acid"),
        ("CCC(C)CC(O)=O", "3-methylvaleric acid"),
        ("OC(=O)\\C=C\\C=C", "(E)-penta-2,4-dienoic acid"),
        ("OC(=O)CC=C", "but-3-enoic acid"),
        ("CC(C)C(O)=O", "isobutyric acid"),
        ("CC(C)[C@@H](C)C(O)=O", "(R)-2,3-dimethylbutyric acid"),
        ("CCCC(O)=O", "butyric acid"),
        ("CCC(C)C(O)=O", "2-methylbutyric acid"),
        ("[H]C(C)=C(C)C(O)=O", "2-methylbut-2-enoic acid"),
        ("[H]C(C=C)=CC(O)=O", "penta-2,4-dienoic acid"),
        ("[H]C(CC)=CC(O)=O", "pent-2-enoic acid"),
        ("CC(C)CC(O)=O", "isovaleric acid"),
        ("[H]\\C(C)=C(/C)C(O)=O", "angelic acid"),
        ("OC(=O)CCC=C", "pent-4-enoic acid"),
        ("[H]C(CC(O)=O)=CC", "pent-3-enoic acid"),
        ("[H]\\C(C)=C\\C(O)=O", "isocrotonic acid"),
        ("CC\\C=C/C(O)=O", "cis-pent-2-enoic acid"),
        ("CC(C)(C)CC(O)=O", "3,3-dimethylbutyric acid"),
        ("CCC(O)=O", "propionic acid"),
        ("CC(C)(C)C(O)=O", "pivalic acid"),
        # (Note: several of the following were flagged as false positives or negatives in the previous attempt)
        ("CCCC(CC=C)C(O)=O", "2-n-Propyl-4-pentenoic acid"),
        ("OC(=O)C(CC)(CC)C", "2-ethyl-2-methyl-butanoic acid"),
        ("OC(=O)CC(C)=C", "Isopropenylacetic acid"),
        ("CCCC(CCC)C(O)=O", "valproic acid"),
        ("OC(=O)C(CC)CC", "2-Ethylbutanoic acid"),
        ("OC(=O)C(C)/C=C/C", "2-Methyl-3-pentenoic acid"),
        ("OC(=O)CC(C)C=C", "3-methyl-4-pentenoic acid"),
        ("OC(=O)C=C", "acrylic acid"),
        ("OC(=O)CC(CC)(C)C", "beta,beta-dimethyl valeric acid"),
        ("CCC\\C(=C/C=C)C(O)=O", "2-Propyl-2,4-pentadienoic acid"),
        ("OC(=O)C(CC=C)(C)C", "2,2-dimethyl-4-pentenoic acid"),
        ("CC(O)=O", "acetic acid"),
        ("C=CCC(C)C(=O)O", "2-Methyl-4-pentenoic acid"),
        ("CCCC(\\C=C/C)C(O)=O", "2-n-Propyl-3-pentenoic acid"),
        ("OC(=O)C(=C(C)C)C", "trimethyl acrylic acid"),
        ("[H]C(O)=O", "formic acid"),
        ("OC(=O)C(=CC(C)C)C", "2,4-Dimethyl-2-pentenoic acid"),
        ("OC(=O)\\C=C\\C(C)C", "4-Methyl-2-pentenoic acid"),
        ("OC(=O)\\C(\\C(C)C)=C\\C", "2-isopropyl trans-crotonic acid"),
        ("OC(=O)/C=C(\\CC)/C", "3-methyl-2Z-pentenoic acid"),
        ("OC(CCCCC)C=C.OC(=O)C", "1-octen-3-ol acetate"),
        ("CCCC(=CCC)C(=O)O", "2-propyl-2-pentenoic acid"),
        # Examples that in previous attempt were false negatives (expected to be accepted)
        ("CCC(CO)C(O)=O", "2-ethylhydracrylic acid"),
        ("[H][C@@]12[C@H](CCN1CC=C2COC(=O)[C@](O)([C@H](C)O)C(C)(C)O)OC(=O)C(\\C)=C/C", "heliosupine"),
        ("O[C@H]([C@@H](CC)C)C(O)=O", "(2R,3R)-2-hydroxy-3-methylpentanoic acid"),
        ("OCCC(O)=O", "3-hydroxypropionic acid"),
    ]
    for smi, name in test_examples:
        res, reason = is_short_chain_fatty_acid(smi)
        print(f"SMILES: {smi}\n  NAME: {name}\n  Classified as: {res}\n  Reason: {reason}\n")