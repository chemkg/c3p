"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: short-chain fatty acid
Definition:
  An aliphatic monocarboxylic acid with a chain length (including the carboxyl carbon) of less than C6.
  In addition, substituents on the backbone (if any) must be hydrocarbon in nature – if any substituent (apart 
  from the carboxyl group itself) contains a heteroatom, the compound is not normally regarded as a short‐chain fatty acid.
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    
    A short-chain fatty acid in this context is defined as an aliphatic monocarboxylic acid
    that contains exactly one carboxyl group, and whose longest continuous carbon chain (starting
    at the carboxyl carbon) has fewer than 6 carbons. Also, any substituent attached to that chain 
    (apart from the carboxyl oxygens on the acid carbon) must be purely hydrocarbon (only C and implicit H).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We require that the molecule be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and is not purely aliphatic"
    
    # Look for exactly one carboxylic acid group: [CX3](=O)[OX2H]
    acid_query = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    acid_matches = mol.GetSubstructMatches(acid_query)
    if len(acid_matches) != 1:
        return False, "Molecule does not have exactly one carboxyl group"
    
    # Identify acid carbon and its acid oxygens.
    acid_match = acid_matches[0]
    acid_carbon = acid_match[0]
    acid_oxygens = set(acid_match[1:])  # these oxygens belong to the acid group.
    
    # Build a graph over carbon atoms obtainable in the molecule.
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # For each carbon atom, record neighboring carbon atoms (via C–C bonds only).
    carbon_graph = {idx: [] for idx in carbon_atoms}
    for idx in carbon_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                nidx = nbr.GetIdx()
                if nidx in carbon_graph:
                    carbon_graph[idx].append(nidx)
    
    if acid_carbon not in carbon_graph:
        return False, "Acid carbon is not present in carbon graph"

    # In an acyclic molecule, the carbon graph is a tree.
    # We enumerate all simple paths (list of atom indices) starting at acid_carbon.
    paths = []
    def dfs_path(current, path):
        # Append a copy of current path as a candidate.
        paths.append(path.copy())
        for nbr in carbon_graph[current]:
            if nbr not in path:  # simple path: no repeats
                path.append(nbr)
                dfs_path(nbr, path)
                path.pop()
    dfs_path(acid_carbon, [acid_carbon])
    
    if not paths:
        return False, "No carbon chain found from acid carbon"
    
    # Find the maximum chain length (number of carbons) among all paths.
    max_len = max(len(p) for p in paths)
    if max_len >= 6:
        return False, f"Longest carbon chain has {max_len} carbons, exceeding the short-chain limit (<6)"
    
    # Filter all candidate backbones that achieve this maximum length.
    candidate_backbones = [set(p) for p in paths if len(p)==max_len]
    
    # Function to check that a branch is purely hydrocarbon.
    # We perform a DFS starting at branch_start (which is a neighbor of a backbone carbon)
    # and ensure that every atom visited is either carbon (atomic num 6) or hydrogen (atomic num 1).
    # We avoid going back into the backbone.
    def branch_is_hydrocarbon(start, backbone_set, visited):
        if start in visited:
            return True
        visited.add(start)
        atom = mol.GetAtomWithIdx(start)
        if atom.GetAtomicNum() not in (6, 1):
            return False
        # Only traverse further if the atom is carbon.
        if atom.GetAtomicNum() == 6:
            for nbr in atom.GetNeighbors():
                nidx = nbr.GetIdx()
                # Do not traverse if part of backbone.
                if nidx in backbone_set:
                    continue
                if nidx in visited:
                    continue
                if nbr.GetAtomicNum() not in (6, 1):
                    return False
                if nbr.GetAtomicNum() == 6:
                    if not branch_is_hydrocarbon(nidx, backbone_set, visited):
                        return False
        return True

    # For each candidate backbone, check whether every substituent off a backbone carbon is hydrocarbon.
    # For the acid carbon, skip the acid oxygens.
    candidate_ok = False
    candidate_reasons = []
    for backbone in candidate_backbones:
        backbone_ok = True
        for idx in backbone:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nidx = nbr.GetIdx()
                # If at acid carbon, skip the oxygens that are part of the acid group.
                if idx == acid_carbon and nidx in acid_oxygens:
                    continue
                # Also skip if the neighbor is part of the chosen backbone.
                if nidx in backbone:
                    continue
                # Now, if the neighbor is not C or H, it makes the branch non-hydrocarbon.
                if nbr.GetAtomicNum() not in (6,1):
                    backbone_ok = False
                    candidate_reasons.append(f"Non-hydrocarbon substituent found on carbon {idx}")
                    break
                # If it is a carbon, then check its entire branch.
                if nbr.GetAtomicNum() == 6:
                    if not branch_is_hydrocarbon(nidx, backbone, set()):
                        backbone_ok = False
                        candidate_reasons.append(f"Branch from carbon {idx} contains a non-hydrocarbon atom")
                        break
            if not backbone_ok:
                break
        if backbone_ok:
            candidate_ok = True
            break

    if candidate_ok:
        return True, f"Molecule is an aliphatic monocarboxylic acid with a carbon chain length of {max_len} (<6) and only hydrocarbon substituents on the backbone"
    else:
        # In case no candidate backbone passed the substituent test, show one of the reasons.
        reason = candidate_reasons[0] if candidate_reasons else "Non-hydrocarbon substituent found attached to the acyl chain"
        return False, reason

# (Optional) Main block to test a few provided SMILES strings.
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
        # False positives
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
        # False negatives (should be rejected)
        ("CCC(CO)C(O)=O", "2-ethylhydracrylic acid"),
        ("[H][C@@]12[C@H](CCN1CC=C2COC(=O)[C@](O)([C@H](C)O)C(C)(C)O)OC(=O)C(\\C)=C/C", "heliosupine"),
        ("O[C@H]([C@@H](CC)C)C(O)=O", "(2R,3R)-2-hydroxy-3-methylpentanoic acid"),
        ("OCCC(O)=O", "3-hydroxypropionic acid"),
    ]
    for smi, name in test_examples:
        res, reason = is_short_chain_fatty_acid(smi)
        print(f"SMILES: {smi}  NAME: {name}\n  Classified as: {res}\n  Reason: {reason}\n")