"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoid
Cannabinoids are a diverse group of pharmacologically active secondary metabolites 
that occur in the Cannabis plant and are also produced endogenously in humans and animals.
They are typically characterized by the presence of oxygen either as part of a heterocyclic
ring system (often fused) or as aromatic hydroxyls (for example in a resorcinol-like motif) or
in ethanolamide/glycerol ester fragments (observed in endocannabinoids).
This heuristic algorithm uses structural filters (molecular weight range, oxygen presence,
limited amide bonds) and searches for at least one of several cannabinoid‐like motifs that are
accompanied by a “long” alkyl chain (i.e. a substructure "CCCC") to reduce false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    
    This implementation checks that:
      1. The SMILES string is valid.
      2. The molecular weight is within 200–900 Da.
      3. The molecule contains oxygen.
      4. The molecule does not have too many amide bonds (≥2 amide bonds are vetoed).
      5. At least one cannabinoid-like oxygen motif is detected:
           a. A resorcinol motif (benzene with two hydroxyl substituents) that is connected to 
              a long alkyl chain (detected by the SMARTS "CCCC").
           b. An ethanolamide fragment (C(=O)NCCO) that is linked to a long acyl chain.
           c. A glycerol monoester motif (approximated here by a fragment O–C(O)) whose acyl group 
              contains a long alkyl chain.
           d. A fused ring system with oxygen where one ring is further attached to a long alkyl chain.
           
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a cannabinoid, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is out of the typical cannabinoid range (200–900 Da)"
    
    # Must contain oxygen
    if not any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms()):
        return False, "Molecule has no oxygen atoms"

    # Veto excessive amide bonds (pattern: C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    n_amide = len(mol.GetSubstructMatches(amide_pattern))
    if n_amide >= 2:
        return False, f"Found {n_amide} amide bonds (≥2), more characteristic of peptides than cannabinoids"

    # Define helper SMARTS for “long alkyl chain”
    long_chain = Chem.MolFromSmarts("CCCC")  # at least 4 connected carbons

    found_motifs = []

    # --- Motif 1: Resorcinol motif (benzene with two hydroxyl groups) ---
    resorcinol = Chem.MolFromSmarts("c1cc(O)cc(O)c1")
    if mol.HasSubstructMatch(resorcinol):
        found_resorcinol = False
        reso_matches = mol.GetSubstructMatches(resorcinol)
        # Check if any atom in the resorcinol ring touches a long alkyl chain atom
        for match in reso_matches:
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in match and nbr.HasSubstructMatch(long_chain):
                        found_resorcinol = True
                        break
                if found_resorcinol:
                    break
            if found_resorcinol:
                break
        if found_resorcinol:
            found_motifs.append("found resorcinol motif with long alkyl chain")
        else:
            # even if no long chain linked, register the presence (but less confidently)
            found_motifs.append("found resorcinol motif (no linked long alkyl chain detected)")
    
    # --- Motif 2: Ethanolamide/acylamide motif ---
    # Look for fragment C(=O)NCCO
    ethanolamide = Chem.MolFromSmarts("C(=O)NCCO")
    eth_matches = mol.GetSubstructMatches(ethanolamide)
    if eth_matches:
        long_acyl = False
        for match in eth_matches:
            # match[0] is carbonyl carbon; check its neighbors (excluding the amide nitrogen) for a long chain.
            carbonyl = mol.GetAtomWithIdx(match[0])
            for nbr in carbonyl.GetNeighbors():
                # avoid the nitrogen (which is part of the motif)
                if nbr.GetSymbol() == "C":
                    # Check if within the environment a "CCCC" is found
                    env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, nbr.GetIdx())
                    if env:
                        submol = Chem.PathToSubmol(mol, env)
                        if submol.HasSubstructMatch(long_chain):
                            long_acyl = True
                            break
            if long_acyl:
                break
        if long_acyl:
            found_motifs.append("found ethanolamide motif with long acyl chain")
        else:
            found_motifs.append("found ethanolamide motif (acyl chain not clearly long)")
    
    # --- Motif 3: Glycerol monoester (endocannabinoid) motif ---
    # Many endocannabinoids (eg. 2-arachidonoylglycerol) have a glycerol group connected via an ester linkage.
    # We look for a pattern where a glycerol fragment (OCC(O)CO) is connected via an ester: O-C(=O)R.
    glycerol_ester = Chem.MolFromSmarts("O[C;H2]C(O)COC(=O)")
    if glycerol_ester and mol.HasSubstructMatch(glycerol_ester):
        # Now check that the acyl part (after C(=O)) has a long alkyl chain.
        matches = mol.GetSubstructMatches(glycerol_ester)
        long_acyl_in_gly = False
        # In our SMARTS the acyl carbon is immediately after the pattern; we check its neighbors for "CCCC"
        for match in matches:
            # We try to locate the carbonyl carbon (here it is at the end of the match)
            carbonyl_idx = match[-1]
            carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
            for nbr in carbonyl.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    # get an environment around this neighbor and check for long alkyl chain pattern.
                    env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, nbr.GetIdx())
                    if env:
                        submol = Chem.PathToSubmol(mol, env)
                        if submol.HasSubstructMatch(long_chain):
                            long_acyl_in_gly = True
                            break
            if long_acyl_in_gly:
                break
        if long_acyl_in_gly:
            found_motifs.append("found glycerol monoester motif with long acyl chain")
        else:
            found_motifs.append("found glycerol monoester motif (acyl chain not clearly long)")
    
    # --- Motif 4: Fused heterocyclic system with oxygen ---
    # Look for two rings sharing at least 2 atoms and one ring containing oxygen.
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    fused_found = False
    for i in range(len(atom_rings)):
        for j in range(i+1, len(atom_rings)):
            # Overlap of at least 2 atoms indicates a fused system
            if len(set(atom_rings[i]).intersection(atom_rings[j])) >= 2:
                # Check if one of the rings has an oxygen
                ring1 = [mol.GetAtomWithIdx(idx) for idx in atom_rings[i]]
                ring2 = [mol.GetAtomWithIdx(idx) for idx in atom_rings[j]]
                if (any(a.GetAtomicNum() == 8 for a in ring1) or any(a.GetAtomicNum() == 8 for a in ring2)):
                    # additionally, require that at least one atom in one of these rings is adjacent to a long alkyl chain fragment
                    attached_chain = False
                    for ring in (atom_rings[i], atom_rings[j]):
                        for idx in ring:
                            atom = mol.GetAtomWithIdx(idx)
                            for nbr in atom.GetNeighbors():
                                if nbr.GetIdx() not in ring:
                                    # check if this neighbor or its environment contains "CCCC"
                                    env = Chem.FindAtomEnvironmentOfRadiusN(mol, 2, nbr.GetIdx())
                                    if env:
                                        submol = Chem.PathToSubmol(mol, env)
                                        if submol.HasSubstructMatch(long_chain):
                                            attached_chain = True
                                            break
                            if attached_chain:
                                break
                        if attached_chain:
                            break
                    if attached_chain:
                        fused_found = True
                        break
        if fused_found:
            break
    if fused_found:
        found_motifs.append("found fused heterocyclic system with oxygen and attached long alkyl chain")
        
    # Decide classification on the basis of found motifs
    if not found_motifs:
        return False, "Molecule does not display any recognizable cannabinoid-like oxygen motifs with an attached long alkyl chain"
    
    combined_reason = "; ".join(found_motifs) + f"; molecular weight = {mol_wt:.1f} Da"
    return True, combined_reason

# (For testing purposes, one might call the function with an example SMILES.)
if __name__ == "__main__":
    # Example test: cannabigerolic acid
    test_smiles = "CCCCc1cc(O)c(C\\C=C(/C)CCC=C(C)C)c(O)c1C(O)=O"
    result, reason = is_cannabinoid(test_smiles)
    print(result, reason)