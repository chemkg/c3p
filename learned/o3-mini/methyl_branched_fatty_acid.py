"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
Definition: any branched-chain fatty acid containing methyl branches only.
The molecule must have a carboxylic acid (or its deprotonated form) group and an alkyl chain (starting at the alpha-carbon)
whose only branches (i.e. substituents off the main chain) are simple methyl groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule (given by a SMILES string) is a methyl-branched fatty acid.
    
    Requirements:
      - Contains a carboxylic acid (or deprotonated acid) group.
      - Has an alkyl chain starting from the alpha-carbon (the carbon directly attached to the carboxyl carbon).
      - Any branch off the main alkyl chain must be a simple methyl group, i.e.:
           * The branch carbon must be sp³ hybridized.
           * It must have exactly one heavy-atom neighbor (its connection to the main chain).
           * It must carry exactly 3 hydrogen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a methyl-branched fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Identify carboxylic acid substructure. We cover -COOH and -COO−.
    acid_smarts = "[CX3](=O)[OX2H1,OX1-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Take the first match; the first atom in the match is the carboxyl carbon.
    carboxyl_idx = acid_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Find the alpha-carbon: a neighboring carbon (atomic num 6) of the carboxyl carbon.
    alpha_atom = None
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            alpha_atom = nbr
            break
    if alpha_atom is None:
        return False, "No alpha carbon found attached to carboxyl group"
    
    # Now, build the main carbon chain (backbone)
    # We want the longest uninterrupted chain of carbons starting at the alpha-carbon.
    # We exclude the carboxyl carbon from this DFS.
    def dfs_longest_path(atom, visited):
        # Append current atom to the path
        current_idx = atom.GetIdx()
        path = visited + [current_idx]
        best_path = path
        # Loop over carbon neighbors only.
        for nbr in atom.GetNeighbors():
            # Only follow carbon atoms
            if nbr.GetAtomicNum() != 6:
                continue
            # Exclude visits to the carboxyl carbon.
            if nbr.GetIdx() == carboxyl_idx:
                continue
            if nbr.GetIdx() in visited:
                continue
            candidate = dfs_longest_path(nbr, path)
            if len(candidate) > len(best_path):
                best_path = candidate
        return best_path
    
    main_chain = dfs_longest_path(alpha_atom, [])
    if len(main_chain) < 1:
        return False, "The alkyl chain starting at the alpha carbon is too short"
    
    # Now, for each carbon in the main chain, check for substituents that are not part of the chain
    # (and also skip the carboxyl carbon).
    branch_found = False
    for atom_idx in main_chain:
        backbone_atom = mol.GetAtomWithIdx(atom_idx)
        # Iterate over neighbors that are carbons.
        for nbr in backbone_atom.GetNeighbors():
            # We only care about carbon substituents.
            if nbr.GetAtomicNum() != 6:
                continue
            # Do not consider atoms that are already in the main chain or the carboxyl carbon.
            if nbr.GetIdx() in main_chain or nbr.GetIdx() == carboxyl_idx:
                continue
            # For this substituent branch, check that it qualifies as a simple methyl group.
            # It must be sp3 hybridized.
            if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                return False, f"Found branch at atom idx {nbr.GetIdx()} that is not sp3 hybridized (likely not a simple methyl group)"
            # It should be connected only to the backbone (only one heavy atom neighbor).
            heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1]
            if len(heavy_neighbors) != 1:
                return False, f"Found branch at atom idx {nbr.GetIdx()} that does not have exactly one heavy-atom neighbor"
            # Check that the branch carbon carries exactly 3 hydrogens.
            if nbr.GetTotalNumHs() != 3:
                return False, f"Found branch at atom idx {nbr.GetIdx()} that does not have exactly 3 hydrogens"
            branch_found = True

    if not branch_found:
        return False, "No methyl branches detected on the fatty acid chain"
    
    return True, "Molecule has a carboxyl group with an alkyl chain that only has simple methyl branches"

# Example usage (testing the provided examples):
if __name__ == '__main__':
    test_smiles_list = [
        "OC(=O)CC(C)=C",  # Isopropenylacetic acid (should pass as methyl-branched)
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methyltriacontanoic acid
        "CCC(C)CCCCCCCCCCCCC(O)=O",  # 14-methylhexadecanoic acid
        "CC(C)CCCCCCCCCC(O)=O",  # isotridecanoic acid
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 24-methylpentacosanoic acid
        "OC(=O)C(CCC)CC",  # alpha-ethyl valeric acid (should fail: branch that is not methyl)
        "CCC(C)CCCCCCCCCCCCCCCCCCC(O)=O",  # 20-methyldocosanoic acid
        "[O-]C(=O)CCC([N+](C)(C)C)C",  # 4-aminovaleric acid betaine (should fail: branch issues)
        "CCC(C)C(O)=O",  # 2-methylbutyric acid
        "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\C(O)=O",  # (2E,6E,10E,14E)-omega-hydroxygeranylgeranic acid
        "CC(C)CCCCCCCCCCCCCCCCCCC(O)=O",  # 20-methylhenicosanoic acid
        "CC[C@@H](C)C(O)=O",  # (R)-2-methylbutyric acid
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 24-methylhexacosanoic acid
        "OC(=O)\\C=C\\C(C)C",  # 4-Methyl-2-pentenoic acid
        "CCC(C)CCCCCCCCC(O)=O",  # 10-methyldodecanoic acid
        "OC(=O)C(C(NC(OC)=O)(C)C)(C)C",  # 3-[(methoxycarbonyl)amino]-2,2,3-trimethylbutanoic acid
        "CC(CCCCCCC/C=C/C(=O)O)C",  # (E)-11-methyldodec-2-enoic acid
        "CCCCCCCCC(C)CC(O)=O",  # 3-methylundecanoic acid
        "OC(CCCCCCCCCCCCCCCC(C)C)=O",  # 17-methyloctadecanoic acid
        "OC(=O)\\C=C\\C(C)(C)C",  # 4,4-dimethyl-2E-pentenoic acid
        "CC(C)=CC(O)=O",  # 3-methylbut-2-enoic acid
        "CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 26-methyloctacosanoic acid
        "OC(=O)CCC/C=C\\CC/C=C\\C/C=C\\CCC(CC)C",  # 16-methyl-octadeca-5Z,9Z,12Z-trienoic acid
        "CC(C)=CCCC(C)=CCCC(C)=CC(O)=O",  # farnesoic acid
        "CCC\\C(=C/CC)C(O)=O",  # 2-n-Propyl-2-pentenoic acid
        "OC(C(CCCCCCCCCCCCCCCC)C)=O",  # 2-methyloctadecanoic acid
        "OC(=O)CC(CC(O)=O)C(C)=C",  # 3-Isopropenylpentanedioic acid
        "O=C(O)/C=C(/CCOC(=O)C)/C",  # Pestalotiopin A
        "OC(=O)\\C(\\C(C)C)=C\\C",  # 2-isopropyl trans-crotonic acid
        "OC(=O)/C=C\\C(C)(C)C",  # 4,4-dimethyl-2Z-pentenoic acid
        "C([C@@H](CCC[C@H](CCC[C@H](CCCC(C)C)C)C)(C)=O)O",  # (2R,6S,10S)-2,6,10,14-pristanic acid
        "OC(=O)/C=C(\\CC)/C",  # 3-methyl-2Z-pentenoic acid
        "CC(C)CCCC(C)CCCC(C)CCC(O)=O",  # 4,8,12-trimethyltridecanoic acid
        "OC(=O)CCC(C)=C",  # 4-methyl-4-pentenoic acid
        "CCCCCCC(C)CCCCCCCCCCC(O)=O",  # 12-methyloctadecanoic acid
        "OC(=O)C(CCC)(C)C",  # alpha,alpha-dimethyl valeric acid
        "OC(=O)C(CC)(CC)C",  # 2-ethyl-2-methyl-butanoic acid
        "CC(C)C(O)=O",  # isobutyric acid
        "CC(C)(C)CC(O)=O",  # 3,3-dimethylbutyric acid
        "OC(=O)C(CCCC(C)C)C",  # 2,6-dimethylheptanoic acid
        "OC(=O)CC(C)C=C",  # 3-methyl-4-pentenoic acid
        "C(O)(=O)CCCCCCCCC(CCC(C)C)C",  # 10,13-dimethyltetradecanoic acid
        "CCCCCC(C)CCCCCC(C)CCCCCCCCCC(C)CC(O)=O",  # 3,13,19-trimethyltricosanoic acid
        "O=C(O)/C=C(\\CC(=O)O)/C",  # Z-3-Methylpent-2-en-1,5-dioic acid
        "CC(CO)CCCCCCCCCCCCCC(O)=O",  # omega-hydroxy-15-methylpalmitic acid
        "OC(=O)CC(CC)(C)C",  # beta,beta-dimethyl valeric acid
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methylnonacosanoic acid
        "O=C(NC(C(CC)C)COC(=O)C)C(CC(=O)O)C",  # 4-((1-acetoxy-3-methylpentan-2-yl)amino)-3-methyl-4-oxobutanoic acid
        "CCC(C)CC(O)=O",  # 3-methylvaleric acid
        "C(C(=CC(=O)[O-])C)CC=C(C)C"  # 3,7-dimethylocta-2,6-dienoate(1-)
    ]
    
    for smi in test_smiles_list:
        result, reason = is_methyl_branched_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")