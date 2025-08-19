
"""
Enrichment analysis of chemical descriptors between class groups where
Chebifier vs Ensemble-9 (rule-like RDKit programs) perform differently.

This version uses **type hints** and the **Typer** CLI, and integrates with 
the c3p data model for easier access to ChEBI class data.

Inputs
------
1) molecules.csv : CSV with columns
   - compound_id : str
   - smiles      : str
   - classes     : str  (semicolon- or pipe-separated list of class names)
   OR a dataset JSON file from the c3p project

2) class_perf.csv : CSV with columns
   - chemical_class : str (exactly matching names used in molecules.csv)
   - f1_left        : float (ensemble-9 F1)   [optional]
   - f1_right       : float (chebifier F1)    [optional]
   - difference     : float (ensemble-9 - chebifier)
   OR can use the pairwise_comparison.csv from c3p-compare output

Outputs
-------
- class_descriptor_table.csv : class-level descriptor table merged with perf labels
- mwut_enrichment.csv        : enrichment stats (U, pval, qval, medians, effect size)
- logreg_coefficients.csv    : standardized logistic regression coefficients

Usage
-----
# Using CSV files:
python analyze_differences.py run \
  --molecules molecules.csv \
  --class-perf class_perf.csv \
  --delta 0.05 \
  --agg median \
  --outdir results/

# Using c3p dataset and comparison output:
python analyze_differences.py explore-differences \
  results/2025/comparison-chebifier/pairwise_comparison.csv \
  --dataset data/dataset.json \
  --delta 0.05 \
  --outdir results/analysis/
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import typer
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from rdkit.Chem.rdMolDescriptors import CalcNumAtomStereoCenters
from scipy.stats import mannwhitneyu
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# Import plotting libraries
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    logging.warning("Matplotlib/Seaborn not available; plotting features will be disabled")

# Import c3p data model
try:
    from c3p.datamodel import Dataset, ChemicalClass, ChemicalStructure
except ImportError:
    Dataset = None
    ChemicalClass = None
    ChemicalStructure = None
    logging.warning("c3p datamodel not available; some features will be limited")

app = typer.Typer(add_completion=False, help="RDKit descriptor enrichment for class-level performance deltas.")
logger = logging.getLogger("rdkit_enrichment")

# ===========================
# SMARTS patterns (typed)
# ===========================
ESTER: Chem.Mol = Chem.MolFromSmarts("C(=O)O")
AMIDE: Chem.Mol = Chem.MolFromSmarts("C(=O)N")
THIOESTER: Chem.Mol = Chem.MolFromSmarts("C(=O)S")
AMINE: Chem.Mol = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]")
QUAT_N: Chem.Mol = Chem.MolFromSmarts("[NX4+]")
PHOSPHATE: Chem.Mol = Chem.MolFromSmarts("P(=O)(O)(O)")
PHOSPHONATE: Chem.Mol = Chem.MolFromSmarts("P(=O)(O)O[C,c]")
ADENINE: Chem.Mol = Chem.MolFromSmarts("n1cnc2ncnc12")
PANT_THIO: Chem.Mol = THIOESTER
EPOXIDE: Chem.Mol = Chem.MolFromSmarts("C1OC1")
ACETAL: Chem.Mol = Chem.MolFromSmarts("[CX4H1](-[OX2])(-[OX2])-[OX2]")

# Enhanced patterns for better analysis
STEROID_NUCLEUS: Chem.Mol = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6](~[#6]~[#6]~[#6]~[#6]~4)~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
GLYCEROL_BACKBONE: Chem.Mol = Chem.MolFromSmarts("[OX2][CX4H2][CX4H]([OX2])[CX4H2][OX2]")
ACYL_TO_GLYCEROL: Chem.Mol = Chem.MolFromSmarts("C(=O)O[CH2,CH][CH2,CH]O")
GLYCOSIDIC_BOND: Chem.Mol = Chem.MolFromSmarts("O[C;R]([O,N])[C;R]")
PHOSPHATE_DIESTER: Chem.Mol = Chem.MolFromSmarts("P(=O)(O[C,c])(O[C,c])")
PHOSPHATE_MONOESTER: Chem.Mol = Chem.MolFromSmarts("P(=O)(O)(O[C,c])")
PYROPHOSPHATE: Chem.Mol = Chem.MolFromSmarts("P(=O)(O)(O)OP(=O)(O)O")
PYRANOSE_RING: Chem.Mol = Chem.MolFromSmarts("C1OCCCC1")
FURANOSE_RING: Chem.Mol = Chem.MolFromSmarts("C1OCCC1")
SIX_MEMBERED_ALIPHATIC: Chem.Mol = Chem.MolFromSmarts("[C;R6]1[C;R6][C;R6][C;R6][C;R6][C;R6]1")

# ===========================
# Dataclasses / config
# ===========================
@dataclass(frozen=True)
class RunConfig:
    molecules: Path
    class_perf: Path
    delta: float = 0.05
    agg: str = "median"  # or "mean"
    outdir: Path = Path("results")
    verbose: bool = False


# ===========================
# Descriptor helpers (typed)
# ===========================

def fused_ring_count(mol: Chem.Mol) -> int:
    ri = mol.GetRingInfo()
    atom_rings = [set(r) for r in ri.AtomRings()]
    fused = 0
    for i in range(len(atom_rings)):
        for j in range(i + 1, len(atom_rings)):
            if len(atom_rings[i].intersection(atom_rings[j])) >= 2:
                fused += 1
    return fused


def longest_aliphatic_chain(mol: Chem.Mol) -> int:
    atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.GetIsAromatic()]
    if not atoms:
        return 0
    adj = GetAdjacencyMatrix(mol)
    longest = 1
    for start in atoms:
        visited: Dict[int, int] = {start: 0}
        queue: List[int] = [start]
        while queue:
            u = queue.pop(0)
            for v, connected in enumerate(adj[u]):
                if connected and v in atoms and v not in visited:
                    visited[v] = visited[u] + 1
                    queue.append(v)
        longest = max(longest, max(visited.values()) + 1)
    return int(longest)


def count_matches(mol: Chem.Mol, patt: Optional[Chem.Mol]) -> int:
    if patt is None:
        return 0
    return len(mol.GetSubstructMatches(patt))


def sugar_like_ring_count(mol: Chem.Mol) -> int:
    ri = mol.GetRingInfo()
    count = 0
    for ring in ri.AtomRings():
        if len(ring) not in (5, 6):
            continue
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if sum(1 for a in ring_atoms if a.GetAtomicNum() in (6, 8)) < len(ring) - 1:
            continue
        substituent_like = 0
        for a in ring_atoms:
            for nbr in a.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    substituent_like += 1
                    break
        if substituent_like >= 3:
            count += 1
    return count


def steroid_like_score(mol: Chem.Mol) -> float:
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    fused = fused_ring_count(mol)
    fr_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    aliph_rings = sum(1 for r in ri.BondRings() if all(not mol.GetBondWithIdx(b).GetIsAromatic() for b in r))
    return 0.5 * fused + 0.3 * aliph_rings + 0.2 * num_rings + 2.0 * fr_sp3


def compute_descriptors(smiles: str) -> Dict[str, float]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"valid": 0.0}
    Chem.SanitizeMol(mol)
    ri = mol.GetRingInfo()
    
    # Calculate bridgehead and spiro atoms
    try:
        num_bridgehead = float(rdMolDescriptors.CalcNumBridgeheadAtoms(mol))
        num_spiro = float(rdMolDescriptors.CalcNumSpiroAtoms(mol))
    except:
        num_bridgehead = 0.0
        num_spiro = 0.0
    
    desc: Dict[str, float] = {
        "valid": 1.0,
        "MolWt": float(Descriptors.MolWt(mol)),
        "NumRings": float(ri.NumRings()),
        "NumAromaticRings": float(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "NumAliphaticRings": float(rdMolDescriptors.CalcNumAliphaticRings(mol)),
        "NumHBD": float(rdMolDescriptors.CalcNumHBD(mol)),
        "NumHBA": float(rdMolDescriptors.CalcNumHBA(mol)),
        "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
        "NumRotatableBonds": float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "FractionCSP3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "NumChiralCenters": float(CalcNumAtomStereoCenters(mol)),
        "LongestAliphaticChain": float(longest_aliphatic_chain(mol)),
        "EsterCount": float(count_matches(mol, ESTER)),
        "AmideCount": float(count_matches(mol, AMIDE)),
        "ThioesterCount": float(count_matches(mol, THIOESTER)),
        "AmineCount": float(count_matches(mol, AMINE)),
        "QuatAmmoniumCount": float(count_matches(mol, QUAT_N)),
        "PhosphateCount": float(count_matches(mol, PHOSPHATE)),
        "PhosphonateLikeCount": float(count_matches(mol, PHOSPHONATE)),
        "EpoxideCount": float(count_matches(mol, EPOXIDE)),
        "AcetalLikeCount": float(count_matches(mol, ACETAL)),
        "SugarLikeRingCount": float(sugar_like_ring_count(mol)),
        "SteroidLikeScore": float(steroid_like_score(mol)),
        "FusedRingPairs": float(fused_ring_count(mol)),
        "HeteroAtomFrac": float(sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6)) / max(1, mol.GetNumAtoms())),
        "CoALikeCue": float(int((count_matches(mol, ADENINE) > 0) and (count_matches(mol, PANT_THIO) > 0) and (count_matches(mol, PHOSPHATE) >= 2))),
        # New enhanced descriptors
        "NumBridgeheadAtoms": num_bridgehead,
        "NumSpiroAtoms": num_spiro,
        "SteroidNucleusCount": float(count_matches(mol, STEROID_NUCLEUS)),
        "GlycerolBackboneCount": float(count_matches(mol, GLYCEROL_BACKBONE)),
        "AcylToGlycerolCount": float(count_matches(mol, ACYL_TO_GLYCEROL)),
        "GlycosidicBondCount": float(count_matches(mol, GLYCOSIDIC_BOND)),
        "PhosphateDiesterCount": float(count_matches(mol, PHOSPHATE_DIESTER)),
        "PhosphateMonoesterCount": float(count_matches(mol, PHOSPHATE_MONOESTER)),
        "PyrophosphateCount": float(count_matches(mol, PYROPHOSPHATE)),
        "PyranoseRingCount": float(count_matches(mol, PYRANOSE_RING)),
        "FuranoseRingCount": float(count_matches(mol, FURANOSE_RING)),
        "SixMemberedAliphaticCount": float(count_matches(mol, SIX_MEMBERED_ALIPHATIC)),
    }
    return desc


# ===========================
# Stats / ML helpers (typed)
# ===========================

def benjamini_hochberg(pvals: Sequence[float]) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = np.asarray(pvals)[order]
    qvals = np.empty(n, dtype=float)
    min_coeff = 1.0
    for i in range(n - 1, -1, -1):
        coeff = n / (i + 1) * ranked[i]
        min_coeff = min(min_coeff, coeff)
        qvals[i] = min_coeff
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(qvals, 1.0)
    return out


def cliffs_delta(x: Sequence[float], y: Sequence[float]) -> float:
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)
    gt = 0
    lt = 0
    for xi in x_arr:
        gt += np.sum(xi > y_arr)
        lt += np.sum(xi < y_arr)
    n = x_arr.size * y_arr.size
    if n == 0:
        return float("nan")
    return float((gt - lt) / n)


# ===========================
# Core pipeline (typed)
# ===========================

def read_tables(molecules: Path, class_perf: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df_mol = pd.read_csv(molecules)
    df_perf = pd.read_csv(class_perf)
    return df_mol, df_perf


def load_dataset(dataset_path: Path) -> Optional[Dict]:
    """Load a c3p dataset from JSON file (returns raw dict for flexibility)."""
    try:
        with open(dataset_path, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        logger.error(f"Failed to load dataset: {e}")
        return None


def dataset_to_molecules_df(dataset: Dict) -> pd.DataFrame:
    """Convert a c3p Dataset dict to a molecules DataFrame format."""
    rows = []
    
    # Create a mapping of SMILES to chemical names
    smiles_to_name = {}
    if "structures" in dataset and dataset["structures"]:
        for struct in dataset["structures"]:
            smiles_to_name[struct["smiles"]] = struct["name"]
    
    # Process each class
    for chem_class in dataset.get("classes", []):
        class_name = chem_class.get("name")
        
        # Handle instances if present (the main field in dataset.json)
        if "instances" in chem_class and chem_class["instances"]:
            for instance in chem_class["instances"]:
                smiles = instance.get("smiles")
                name = instance.get("name", smiles_to_name.get(smiles, f"compound_{len(rows)}"))
                if smiles:  # Only add if SMILES is present
                    rows.append({
                        "compound_id": name,
                        "smiles": smiles,
                        "classes": class_name
                    })
        
        # Handle all_positive_examples if present
        elif "all_positive_examples" in chem_class and chem_class["all_positive_examples"]:
            for smiles in chem_class["all_positive_examples"]:
                name = smiles_to_name.get(smiles, f"compound_{len(rows)}")
                rows.append({
                    "compound_id": name,
                    "smiles": smiles,
                    "classes": class_name
                })
    
    return pd.DataFrame(rows)


def detect_class_column(df_mol: pd.DataFrame) -> str:
    for c in ("classes", "class", "chemical_class", "labels"):
        if c in df_mol.columns:
            return c
    raise ValueError("molecules.csv must contain a 'classes' (or 'class'/'chemical_class'/'labels') column")


def build_descriptor_table(df_mol: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, float]] = []
    for idx, row in df_mol.iterrows():
        smi = str(row["smiles"])  # type: ignore[index]
        d = compute_descriptors(smi)
        d.update({"compound_id": row.get("compound_id", idx), "smiles": smi})
        rows.append(d)
    return pd.DataFrame(rows)


def expand_classes(df_mol: pd.DataFrame, class_col: str) -> pd.DataFrame:
    exp_rows: List[Dict[str, str]] = []
    series = df_mol[class_col].astype(str)
    sep = ";" if series.str.contains(";").any() else "|"
    for idx, row in df_mol.iterrows():
        classes = str(row[class_col]).split(sep)
        for cl in [c.strip() for c in classes if c.strip()]:
            exp_rows.append({"compound_id": row.get("compound_id", idx), "class": cl})
    return pd.DataFrame(exp_rows)


def aggregate_class_descriptors(df_expanded: pd.DataFrame, df_desc: pd.DataFrame, agg: str) -> pd.DataFrame:
    df_desc2 = df_desc.copy()
    # Convert index to series for fillna
    df_desc2["compound_id"] = df_desc2["compound_id"].fillna(pd.Series(df_desc2.index, index=df_desc2.index))
    df_cls_mol = df_expanded.merge(df_desc2, on="compound_id", how="left")

    agg_func: Dict[str, str] = {}
    for c in df_desc.columns:
        if c not in ("compound_id", "smiles"):
            agg_func[c] = "mean" if c == "valid" else agg
    df_cls = df_cls_mol.groupby("class").agg(agg_func).reset_index()
    return df_cls


def label_groups(df_cls: pd.DataFrame, df_perf: pd.DataFrame, delta: float) -> pd.DataFrame:
    df_perf2 = df_perf.rename(columns={"chemical_class": "class"})
    df = df_cls.merge(df_perf2[["class", "difference"]], on="class", how="inner")
    df["group"] = np.where(df["difference"] > delta, "ensemble_favored",
                            np.where(df["difference"] < -delta, "chebifier_favored", "neutral"))
    return df


def run_mwut(df_cls: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    df_fav = df_cls[df_cls["group"].isin(["ensemble_favored", "chebifier_favored"])].copy()
    desc_cols = [c for c in df_cls.columns if c not in ("class", "difference", "group")]

    stats_rows: List[Dict[str, float]] = []
    for col in desc_cols:
        x = df_fav[df_fav["group"] == "ensemble_favored"][col].dropna().values
        y = df_fav[df_fav["group"] == "chebifier_favored"][col].dropna().values
        if len(x) == 0 or len(y) == 0:
            continue
        try:
            u, p = mannwhitneyu(x, y, alternative="two-sided")
        except ValueError:
            u, p = np.nan, 1.0
        eff = cliffs_delta(x, y)
        stats_rows.append({
            "descriptor": col,
            "median_ensemble": float(np.median(x)) if len(x) else float("nan"),
            "median_chebifier": float(np.median(y)) if len(y) else float("nan"),
            "U_stat": float(u) if not np.isnan(u) else float("nan"),
            "p_value": float(p),
            "cliffs_delta": float(eff),
            "n_ensemble": int(len(x)),
            "n_chebifier": int(len(y)),
        })

    df_stats = pd.DataFrame(stats_rows)
    if not df_stats.empty:
        df_stats["q_value"] = benjamini_hochberg(df_stats["p_value"].values)
        df_stats.sort_values(["q_value", "p_value"], inplace=True)
        df_stats.to_csv(out_dir / "mwut_enrichment.csv", index=False)
    return df_stats


def create_enrichment_plots(df_stats: pd.DataFrame, out_dir: Path) -> None:
    """Create visualization plots for enrichment analysis."""
    if not PLOTTING_AVAILABLE:
        logger.warning("Plotting libraries not available, skipping plots")
        return
    
    if df_stats.empty:
        logger.warning("No enrichment data to plot")
        return
    
    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = (12, 8)
    
    # Create a figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Volcano plot (effect size vs significance)
    ax1 = axes[0, 0]
    df_plot = df_stats.copy()
    df_plot['-log10(q)'] = -np.log10(df_plot['q_value'].clip(lower=1e-10))
    
    # Color points by significance
    colors = ['red' if q < 0.05 else 'gray' for q in df_plot['q_value']]
    
    ax1.scatter(df_plot['cliffs_delta'], df_plot['-log10(q)'], 
                c=colors, alpha=0.6, s=50)
    ax1.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='q=0.05')
    ax1.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    ax1.axvline(x=0.33, color='blue', linestyle='--', alpha=0.3, label='|δ|=0.33 (medium)')
    ax1.axvline(x=-0.33, color='blue', linestyle='--', alpha=0.3)
    
    # Label significant points
    sig_points = df_plot[df_plot['q_value'] < 0.05]
    for _, row in sig_points.iterrows():
        if abs(row['cliffs_delta']) > 0.2:  # Only label if effect size is meaningful
            ax1.annotate(row['descriptor'], 
                        (row['cliffs_delta'], row['-log10(q)']),
                        fontsize=8, alpha=0.7)
    
    ax1.set_xlabel("Cliff's Delta (Effect Size)\n← Favors Chebifier | Favors Ensemble →")
    ax1.set_ylabel("-log10(q-value)")
    ax1.set_title("Volcano Plot: Descriptor Enrichment")
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # 2. Top enriched descriptors bar plot
    ax2 = axes[0, 1]
    top_n = 15
    df_sig = df_stats[df_stats['q_value'] < 0.1].sort_values('q_value').head(top_n)
    
    if not df_sig.empty:
        y_pos = np.arange(len(df_sig))
        colors_bar = ['#2E86AB' if d > 0 else '#A23B72' for d in df_sig['cliffs_delta']]
        
        ax2.barh(y_pos, df_sig['cliffs_delta'], color=colors_bar, alpha=0.7)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(df_sig['descriptor'], fontsize=9)
        ax2.set_xlabel("Cliff's Delta (Effect Size)")
        ax2.set_title(f"Top {len(df_sig)} Enriched Descriptors (q < 0.1)")
        ax2.axvline(x=0, color='black', linestyle='-', alpha=0.5)
        ax2.grid(True, alpha=0.3, axis='x')
    else:
        ax2.text(0.5, 0.5, 'No significantly enriched descriptors', 
                ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title("Top Enriched Descriptors")
    
    # 3. Median values comparison (top descriptors)
    ax3 = axes[1, 0]
    df_top = df_stats.sort_values('q_value').head(10)
    
    if not df_top.empty:
        x = np.arange(len(df_top))
        width = 0.35
        
        bars1 = ax3.bar(x - width/2, df_top['median_ensemble'], width, 
                       label='C3P-favored classes', color='#2E86AB', alpha=0.7)
        bars2 = ax3.bar(x + width/2, df_top['median_chebifier'], width,
                       label='Chebifier-favored classes', color='#A23B72', alpha=0.7)
        
        ax3.set_xlabel('Descriptor')
        ax3.set_ylabel('Median Value')
        ax3.set_title('Median Descriptor Values\n(Classes where C3P excels vs Classes where Chebifier excels)')
        ax3.set_xticks(x)
        ax3.set_xticklabels(df_top['descriptor'], rotation=45, ha='right', fontsize=9)
        ax3.legend()
        ax3.grid(True, alpha=0.3, axis='y')
    else:
        ax3.text(0.5, 0.5, 'No data to display', 
                ha='center', va='center', transform=ax3.transAxes)
    
    # 4. Q-value distribution
    ax4 = axes[1, 1]
    ax4.hist(df_stats['q_value'], bins=30, color='steelblue', alpha=0.7, edgecolor='black')
    ax4.axvline(x=0.05, color='red', linestyle='--', alpha=0.7, label='q=0.05')
    ax4.axvline(x=0.1, color='orange', linestyle='--', alpha=0.7, label='q=0.1')
    ax4.set_xlabel('Q-value (FDR-adjusted)')
    ax4.set_ylabel('Count')
    ax4.set_title('Distribution of Q-values')
    ax4.legend()
    ax4.grid(True, alpha=0.3, axis='y')
    
    # Add info text
    n_sig_005 = len(df_stats[df_stats['q_value'] < 0.05])
    n_sig_01 = len(df_stats[df_stats['q_value'] < 0.1])
    ax4.text(0.95, 0.95, f'q<0.05: {n_sig_005}\nq<0.1: {n_sig_01}', 
            transform=ax4.transAxes, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle('Molecular Descriptor Enrichment Analysis\nEnsemble vs Chebifier Performance', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save the plot
    plot_path = out_dir / "enrichment_analysis.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    logger.info(f"Saved enrichment plots to {plot_path}")
    plt.close()
    
    # Create a separate heatmap for top descriptors
    create_descriptor_heatmap(df_stats, df_cls=None, out_dir=out_dir)


def create_descriptor_heatmap(df_stats: pd.DataFrame, df_cls: Optional[pd.DataFrame], out_dir: Path) -> None:
    """Create a heatmap showing descriptor patterns."""
    if not PLOTTING_AVAILABLE:
        return
    
    # Select top descriptors by significance
    top_desc = df_stats.sort_values('q_value').head(20)['descriptor'].tolist()
    
    if len(top_desc) < 2:
        logger.warning("Not enough significant descriptors for heatmap")
        return
    
    # If we have the class descriptor table, use it for the heatmap
    cls_desc_file = out_dir / "class_descriptor_table.csv"
    if cls_desc_file.exists():
        df_cls = pd.read_csv(cls_desc_file)
        
        # Filter to ensemble/chebifier favored classes
        df_fav = df_cls[df_cls['group'].isin(['ensemble_favored', 'chebifier_favored'])]
        
        if not df_fav.empty and len(top_desc) > 0:
            # Prepare data for heatmap
            desc_cols = [col for col in top_desc if col in df_fav.columns]
            
            if desc_cols:
                plt.figure(figsize=(12, 8))
                
                # Normalize descriptor values for better visualization
                from sklearn.preprocessing import StandardScaler
                scaler = StandardScaler()
                df_norm = df_fav[desc_cols].fillna(0)
                df_scaled = pd.DataFrame(
                    scaler.fit_transform(df_norm),
                    columns=desc_cols,
                    index=df_fav.index
                )
                
                # Add class names and groups
                df_scaled['class'] = df_fav['class'].values
                df_scaled['group'] = df_fav['group'].values
                
                # Sort by group
                df_scaled = df_scaled.sort_values('group')
                
                # Create heatmap
                sns.heatmap(df_scaled[desc_cols].T, 
                           cmap='RdBu_r', center=0,
                           cbar_kws={'label': 'Standardized Value'},
                           xticklabels=False,
                           yticklabels=desc_cols)
                
                # Add group labels
                group_positions = []
                current_group = None
                for i, group in enumerate(df_scaled['group']):
                    if group != current_group:
                        group_positions.append((i, group))
                        current_group = group
                
                for pos, group in group_positions:
                    color = '#2E86AB' if 'ensemble' in group else '#A23B72'
                    plt.axvline(x=pos, color=color, linewidth=2, alpha=0.5)
                
                plt.title('Descriptor Patterns: Ensemble vs Chebifier Favored Classes')
                plt.xlabel('Chemical Classes (sorted by group)')
                plt.ylabel('Molecular Descriptors')
                plt.tight_layout()
                
                heatmap_path = out_dir / "descriptor_heatmap.png"
                plt.savefig(heatmap_path, dpi=150, bbox_inches='tight')
                logger.info(f"Saved descriptor heatmap to {heatmap_path}")
                plt.close()


def generate_interpretation_summary(df_stats: pd.DataFrame, regression_results: pd.DataFrame, out_dir: Path) -> None:
    """
    Generate a human-readable interpretation summary of the analysis.
    """
    summary_lines = []
    summary_lines.append("=" * 80)
    summary_lines.append("MOLECULAR DESCRIPTOR ENRICHMENT ANALYSIS SUMMARY")
    summary_lines.append("Ensemble-9 vs Chebifier Performance Differences")
    summary_lines.append("=" * 80)
    summary_lines.append("")
    
    # Key findings based on FDR-corrected significance
    summary_lines.append("KEY FINDINGS (FDR-corrected):")
    summary_lines.append("-" * 40)
    
    if not df_stats.empty:
        sig_005 = df_stats[df_stats['q_value'] < 0.05].sort_values('q_value')
        sig_010 = df_stats[df_stats['q_value'] < 0.10].sort_values('q_value')
        
        if not sig_005.empty:
            summary_lines.append("\n✓ SIGNIFICANT DESCRIPTORS (q < 0.05):")
            for _, row in sig_005.iterrows():
                effect_size = abs(row['cliffs_delta'])
                if effect_size < 0.147:
                    effect = "negligible"
                elif effect_size < 0.33:
                    effect = "small"
                elif effect_size < 0.474:
                    effect = "medium"
                else:
                    effect = "large"
                
                direction = "CHEBIFIER" if row['cliffs_delta'] < 0 else "ENSEMBLE"
                summary_lines.append(f"  • {row['descriptor']}: Higher in {direction}")
                summary_lines.append(f"    - Median values: Ensemble={row['median_ensemble']:.2f}, Chebifier={row['median_chebifier']:.2f}")
                summary_lines.append(f"    - Effect size: Cliff's δ={row['cliffs_delta']:.3f} ({effect})")
                summary_lines.append(f"    - Q-value: {row['q_value']:.4f}")
                summary_lines.append("")
        
        if len(sig_010) > len(sig_005):
            summary_lines.append("\n⚠ BORDERLINE SIGNIFICANT (0.05 < q < 0.10):")
            borderline = sig_010[~sig_010.index.isin(sig_005.index)]
            for _, row in borderline.iterrows():
                direction = "chebifier" if row['cliffs_delta'] < 0 else "ensemble"
                summary_lines.append(f"  • {row['descriptor']}: Trends higher for {direction} (q={row['q_value']:.3f})")
    
    # Interpretation based on descriptor patterns
    summary_lines.append("\n\nINTERPRETATION:")
    summary_lines.append("-" * 40)
    
    # Check for ring-related patterns
    ring_descriptors = ['NumRings', 'NumAromaticRings', 'NumAliphaticRings', 'FusedRingPairs', 
                       'NumBridgeheadAtoms', 'NumSpiroAtoms', 'SteroidNucleusCount', 'SixMemberedAliphaticCount']
    ring_sigs = df_stats[df_stats['descriptor'].isin(ring_descriptors) & (df_stats['q_value'] < 0.10)]
    
    if not ring_sigs.empty:
        summary_lines.append("\n📍 RING COMPLEXITY:")
        for _, row in ring_sigs.iterrows():
            if row['cliffs_delta'] < 0:
                summary_lines.append(f"  - {row['descriptor']}: Chebifier performs better with more complex ring systems")
            else:
                summary_lines.append(f"  - {row['descriptor']}: Ensemble performs better with simpler structures")
    
    # Check for lipid-related patterns
    lipid_descriptors = ['EsterCount', 'GlycerolBackboneCount', 'AcylToGlycerolCount', 
                        'PhosphateCount', 'PhosphateDiesterCount', 'PhosphateMonoesterCount']
    lipid_sigs = df_stats[df_stats['descriptor'].isin(lipid_descriptors) & (df_stats['q_value'] < 0.15)]
    
    if not lipid_sigs.empty:
        summary_lines.append("\n📍 LIPID FEATURES:")
        for _, row in lipid_sigs.iterrows():
            if row['cliffs_delta'] > 0:
                summary_lines.append(f"  - {row['descriptor']}: Ensemble excels with lipid-like features (q={row['q_value']:.3f})")
    
    # Check for carbohydrate patterns
    sugar_descriptors = ['GlycosidicBondCount', 'PyranoseRingCount', 'FuranoseRingCount', 
                        'SugarLikeRingCount', 'AcetalLikeCount']
    sugar_sigs = df_stats[df_stats['descriptor'].isin(sugar_descriptors) & (df_stats['q_value'] < 0.15)]
    
    if not sugar_sigs.empty:
        summary_lines.append("\n📍 CARBOHYDRATE FEATURES:")
        for _, row in sugar_sigs.iterrows():
            summary_lines.append(f"  - {row['descriptor']}: q={row['q_value']:.3f}, δ={row['cliffs_delta']:.3f}")
    elif df_stats['descriptor'].isin(sugar_descriptors).any():
        summary_lines.append("\n📍 CARBOHYDRATE FEATURES: No significant enrichment detected")
        summary_lines.append("  - Current SMARTS patterns may be too conservative")
    
    # Add regression insights if available
    if not regression_results.empty:
        summary_lines.append("\n\n📊 CONTINUOUS REGRESSION INSIGHTS:")
        summary_lines.append("-" * 40)
        top_features = regression_results[regression_results['abs_coefficient'] > 0.01].head(5)
        for _, row in top_features.iterrows():
            direction = "ensemble" if row['coefficient'] > 0 else "chebifier"
            summary_lines.append(f"  • {row['descriptor']}: Strong predictor for {direction} (coef={row['coefficient']:.3f})")
    
    # Statistical summary
    summary_lines.append("\n\nSTATISTICAL SUMMARY:")
    summary_lines.append("-" * 40)
    summary_lines.append(f"  Total descriptors tested: {len(df_stats)}")
    summary_lines.append(f"  Significant after FDR (q<0.05): {len(df_stats[df_stats['q_value'] < 0.05])}")
    summary_lines.append(f"  Borderline significant (q<0.10): {len(df_stats[df_stats['q_value'] < 0.10])}")
    
    # Recommendations
    summary_lines.append("\n\nRECOMMENDATIONS:")
    summary_lines.append("-" * 40)
    summary_lines.append("  1. Consider using ensemble methods for lipid-rich and ester-containing molecules")
    summary_lines.append("  2. Consider using chebifier/DL for polycyclic and ring-complex structures")
    summary_lines.append("  3. Refine carbohydrate detection patterns for better discrimination")
    summary_lines.append("  4. Consider weighted analysis by class size to reduce noise from small classes")
    
    # Write summary to file
    summary_text = "\n".join(summary_lines)
    summary_path = out_dir / "analysis_interpretation.txt"
    with open(summary_path, 'w') as f:
        f.write(summary_text)
    
    logger.info(f"Saved interpretation summary to {summary_path}")
    
    # Also print key findings to console
    print("\n" + "="*60)
    print("ANALYSIS HIGHLIGHTS:")
    print("="*60)
    if not sig_005.empty:
        print(f"✓ Found {len(sig_005)} significant descriptor(s) after FDR correction")
        for _, row in sig_005.head(3).iterrows():
            direction = "chebifier" if row['cliffs_delta'] < 0 else "ensemble"
            print(f"  • {row['descriptor']}: favors {direction} (q={row['q_value']:.3f})")
    else:
        print("⚠ No descriptors survived FDR correction at q<0.05")
        if not sig_010.empty:
            print(f"  But {len(sig_010)} showed trends at q<0.10")


def run_continuous_regression(df_cls: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """
    Run continuous regression analysis on performance differences.
    Uses elastic net and returns feature importance.
    """
    from sklearn.linear_model import ElasticNetCV
    from sklearn.preprocessing import RobustScaler
    
    # Filter to classes with valid performance data
    df_analysis = df_cls.dropna(subset=['difference'])
    
    if len(df_analysis) < 10:
        logger.warning("Not enough data for regression analysis")
        return pd.DataFrame()
    
    desc_cols = [c for c in df_cls.columns if c not in ("class", "difference", "group", "chemical_class")]
    
    # Prepare data
    X = df_analysis[desc_cols].fillna(0).values
    y = df_analysis['difference'].values  # Continuous target
    
    # Use robust scaling for outliers
    scaler = RobustScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Elastic Net with cross-validation
    model = ElasticNetCV(
        l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99],
        cv=5,
        max_iter=2000,
        random_state=42
    )
    
    try:
        model.fit(X_scaled, y)
        
        # Get coefficients
        coefs = model.coef_
        
        # Create results dataframe
        results = pd.DataFrame({
            'descriptor': desc_cols,
            'coefficient': coefs,
            'abs_coefficient': np.abs(coefs)
        }).sort_values('abs_coefficient', ascending=False)
        
        # Add interpretation
        results['direction'] = results['coefficient'].apply(
            lambda x: 'favors_ensemble' if x > 0 else 'favors_chebifier' if x < 0 else 'neutral'
        )
        
        # Save results
        results.to_csv(out_dir / "continuous_regression_results.csv", index=False)
        logger.info(f"Saved continuous regression results (R²={model.score(X_scaled, y):.3f})")
        
        return results
        
    except Exception as e:
        logger.error(f"Continuous regression failed: {e}")
        return pd.DataFrame()


def run_logreg(df_cls: pd.DataFrame, out_dir: Path) -> None:
    df_fav = df_cls[df_cls["group"].isin(["ensemble_favored", "chebifier_favored"])].copy()
    desc_cols = [c for c in df_cls.columns if c not in ("class", "difference", "group")]
    y = (df_fav["group"] == "ensemble_favored").astype(int).values
    X = df_fav[desc_cols].values
    pipe = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler(with_mean=True, with_std=True)),
        ("clf", LogisticRegression(max_iter=500, penalty="l2", C=1.0)),
    ])
    try:
        pipe.fit(X, y)
        clf = pipe.named_steps["clf"]
        coefs = clf.coef_.ravel()
        coef_df = pd.DataFrame({"descriptor": desc_cols, "coef": coefs}).sort_values("coef", ascending=False)
        coef_df.to_csv(out_dir / "logreg_coefficients.csv", index=False)
    except Exception as e:  # noqa: BLE001
        (out_dir / "logreg_coefficients.csv").write_text(f"# Logistic regression failed: {e}", encoding="utf-8")


# ===========================
# Typer command
# ===========================
@app.command()
def run(
    molecules: Path = typer.Option(..., exists=True, readable=True, help="Path to molecules.csv"),
    class_perf: Path = typer.Option(..., exists=True, readable=True, help="Path to class_perf.csv"),
    delta: float = typer.Option(0.05, help="|difference| threshold for favored groups"),
    agg: str = typer.Option("median", help="Class-level aggregation", case_sensitive=False),
    outdir: Path = typer.Option(Path("results"), help="Output directory"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable debug logging"),
) -> None:
    """Run the enrichment analysis with typed CLI."""
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="[%(levelname)s] %(message)s")
    cfg = RunConfig(molecules=molecules, class_perf=class_perf, delta=delta, agg=agg.lower(), outdir=outdir, verbose=verbose)

    cfg.outdir.mkdir(parents=True, exist_ok=True)

    df_mol, df_perf = read_tables(cfg.molecules, cfg.class_perf)
    class_col = detect_class_column(df_mol)

    df_desc = build_descriptor_table(df_mol)
    df_expanded = expand_classes(df_mol, class_col)
    df_cls = aggregate_class_descriptors(df_expanded, df_desc, cfg.agg)
    df_labeled = label_groups(df_cls, df_perf, cfg.delta)

    # Save class-level descriptor table
    (cfg.outdir / "class_descriptor_table.csv").write_text("", encoding="utf-8")
    df_labeled.to_csv(cfg.outdir / "class_descriptor_table.csv", index=False)

    # Stats + ML
    _ = run_mwut(df_labeled, cfg.outdir)
    run_logreg(df_labeled, cfg.outdir)

    typer.echo(f"Saved results to: {cfg.outdir.resolve()}")


@app.command()
def analyze_individual_method(
    comparison_file: Path = typer.Argument(..., help="Path to pairwise_comparison.csv from c3p-compare"),
    method: str = typer.Argument(..., help="Method to analyze: 'c3p' or 'chebifier'"),
    dataset: Optional[Path] = typer.Option(None, "--dataset", "-d", help="Path to c3p dataset JSON file"),
    molecules: Optional[Path] = typer.Option(None, "--molecules", "-m", help="Alternative: Path to molecules CSV"),
    f1_threshold: float = typer.Option(0.7, "--threshold", "-t", help="F1 threshold for 'good' performance"),
    agg: str = typer.Option("median", help="Class-level aggregation", case_sensitive=False),
    outdir: Path = typer.Option(Path("results/analysis"), help="Output directory"),
    max_molecules_per_class: Optional[int] = typer.Option(None, "--max-molecules", help="Max molecules per class"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable debug logging"),
) -> None:
    """
    Analyze what molecular features correlate with good/poor performance for a single method.
    
    This identifies molecular descriptors that distinguish classes where a method performs
    well (F1 > threshold) versus poorly (F1 < threshold).
    
    Examples:
        # Analyze c3p (ensemble) performance
        c3p-analyze analyze-individual-method \\
            results/2025/comparison-chebifier/pairwise_comparison.csv c3p
        
        # Analyze chebifier performance
        c3p-analyze analyze-individual-method \\
            results/2025/comparison-chebifier/pairwise_comparison.csv chebifier
    """
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="[%(levelname)s] %(message)s")
    
    # Validate method name
    method_lower = method.lower()
    if method_lower not in ['c3p', 'chebifier', 'ensemble', 'ensemble-9']:
        typer.echo(f"Error: Method must be 'c3p' or 'chebifier', got '{method}'")
        raise typer.Exit(1)
    
    # Normalize method names
    if method_lower in ['ensemble', 'ensemble-9', 'c3p']:
        method_lower = 'c3p'
        method_display = 'C3P (Ensemble-9)'
        f1_column = 'f1_left'  # ensemble is on the left in comparison files
    else:
        method_display = 'Chebifier'
        f1_column = 'f1_right'  # chebifier is on the right

    print(f"Analyzing {method_display} performance; f1_column: {f1_column}")
    
    outdir = outdir / f"individual_{method_lower}"
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Read the comparison file
    df_perf = pd.read_csv(comparison_file)
    
    # Handle different possible column names
    if "chemical_class" not in df_perf.columns and "class" in df_perf.columns:
        df_perf = df_perf.rename(columns={"class": "chemical_class"})
    
    # Get molecules data (reuse existing logic)
    df_mol = None
    if dataset and dataset.exists():
        logger.info(f"Loading dataset from {dataset}")
        ds = load_dataset(dataset)
        if ds:
            df_mol = dataset_to_molecules_df(ds)
            logger.info(f"Loaded {len(df_mol)} molecule-class pairs from dataset")
    elif molecules and molecules.exists():
        df_mol = pd.read_csv(molecules)
    else:
        # Try to find dataset files automatically
        search_paths = [
            Path("results/2025/benchmark/dataset.json"),
            Path("results/2024/inputs/dataset.json"),
        ]
        for path in search_paths:
            if path.exists():
                logger.info(f"Found dataset at {path}")
                ds = load_dataset(path)
                if ds:
                    df_mol = dataset_to_molecules_df(ds)
                    break
    
    if df_mol is None or df_mol.empty:
        logger.error("No molecule data found. Please provide --dataset or --molecules")
        raise typer.Exit(1)
    
    # Detect class column
    class_col = detect_class_column(df_mol)
    
    # Filter to classes in performance data
    perf_classes = set(df_perf["chemical_class"].unique())
    df_mol = df_mol[df_mol[class_col].isin(perf_classes)]
    
    # Limit molecules per class if requested
    if max_molecules_per_class:
        logger.info(f"Limiting to {max_molecules_per_class} molecules per class")
        sampled_rows = []
        for cls in perf_classes:
            cls_rows = df_mol[df_mol[class_col] == cls]
            if len(cls_rows) > max_molecules_per_class:
                cls_rows = cls_rows.sample(n=max_molecules_per_class, random_state=42)
            sampled_rows.append(cls_rows)
        if sampled_rows:
            df_mol = pd.concat(sampled_rows, ignore_index=True)
        logger.info(f"Reduced to {len(df_mol)} molecules total")
    
    # Build descriptors
    logger.info("Computing molecular descriptors...")
    df_desc = build_descriptor_table(df_mol)
    
    # Expand classes and aggregate
    df_expanded = expand_classes(df_mol, class_col)
    df_cls = aggregate_class_descriptors(df_expanded, df_desc, agg.lower())
    
    # Add performance data for the specific method
    df_perf_method = df_perf[["chemical_class", f1_column]].copy()
    df_perf_method.columns = ["class", "f1_score"]
    df_cls = df_cls.merge(df_perf_method, on="class", how="inner")
    
    # Label classes as good/poor performance
    df_cls["performance_group"] = pd.cut(
        df_cls["f1_score"],
        bins=[0, 0.5, f1_threshold, 1.0],
        labels=["poor", "moderate", "good"]
    )
    
    # Save the labeled data
    df_cls.to_csv(outdir / f"{method_lower}_class_descriptors.csv", index=False)
    logger.info(f"Saved class descriptors to {outdir / f'{method_lower}_class_descriptors.csv'}")
    
    # Run MWU tests comparing good vs poor performance
    logger.info(f"Running statistical tests for {method_display}...")
    df_stats = run_method_performance_analysis(df_cls, method_lower, outdir)
    
    # Create visualizations
    if PLOTTING_AVAILABLE:
        create_method_performance_plots(df_cls, df_stats, method_display, outdir)
    
    # Generate interpretation
    generate_method_interpretation(df_cls, df_stats, method_display, f1_threshold, outdir)
    
    typer.echo(f"\n✓ Analysis complete for {method_display}. Results saved to: {outdir.resolve()}")


def run_method_performance_analysis(df_cls: pd.DataFrame, method: str, out_dir: Path) -> pd.DataFrame:
    """
    Run Mann-Whitney U tests comparing molecular descriptors between good and poor performance classes
    for a single method.
    """
    # Filter to good and poor performance groups
    df_analysis = df_cls[df_cls["performance_group"].isin(["good", "poor"])].copy()
    
    if df_analysis.empty:
        logger.warning("No classes in good/poor performance groups")
        return pd.DataFrame()
    
    desc_cols = [c for c in df_cls.columns if c not in ("class", "f1_score", "performance_group")]
    
    stats_rows: List[Dict[str, float]] = []
    for col in desc_cols:
        good = df_analysis[df_analysis["performance_group"] == "good"][col].dropna().values
        poor = df_analysis[df_analysis["performance_group"] == "poor"][col].dropna().values
        
        if len(good) == 0 or len(poor) == 0:
            continue
        
        try:
            u, p = mannwhitneyu(good, poor, alternative="two-sided")
        except ValueError:
            u, p = np.nan, 1.0
        
        eff = cliffs_delta(good, poor)
        
        stats_rows.append({
            "descriptor": col,
            "median_good_perf": float(np.median(good)) if len(good) else float("nan"),
            "median_poor_perf": float(np.median(poor)) if len(poor) else float("nan"),
            "U_stat": float(u) if not np.isnan(u) else float("nan"),
            "p_value": float(p),
            "cliffs_delta": float(eff),
            "n_good": int(len(good)),
            "n_poor": int(len(poor)),
        })
    
    df_stats = pd.DataFrame(stats_rows)
    if not df_stats.empty:
        df_stats["q_value"] = benjamini_hochberg(df_stats["p_value"].values)
        df_stats.sort_values(["q_value", "p_value"], inplace=True)
        df_stats.to_csv(out_dir / f"{method}_performance_enrichment.csv", index=False)
        logger.info(f"Saved performance analysis to {out_dir / f'{method}_performance_enrichment.csv'}")
    
    return df_stats


def create_method_performance_plots(df_cls: pd.DataFrame, df_stats: pd.DataFrame, 
                                   method_name: str, out_dir: Path) -> None:
    """Create visualizations for single method performance analysis."""
    if not PLOTTING_AVAILABLE:
        return
    
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Performance distribution
    ax1 = axes[0, 0]
    ax1.hist(df_cls["f1_score"], bins=30, color='steelblue', alpha=0.7, edgecolor='black')
    ax1.axvline(x=0.5, color='red', linestyle='--', alpha=0.7, label='Poor threshold')
    ax1.axvline(x=0.7, color='green', linestyle='--', alpha=0.7, label='Good threshold')
    ax1.set_xlabel('F1 Score')
    ax1.set_ylabel('Number of Classes')
    ax1.set_title(f'{method_name} Performance Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add statistics
    mean_f1 = df_cls["f1_score"].mean()
    median_f1 = df_cls["f1_score"].median()
    ax1.text(0.02, 0.98, f'Mean: {mean_f1:.3f}\nMedian: {median_f1:.3f}',
            transform=ax1.transAxes, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 2. Volcano plot for descriptors
    ax2 = axes[0, 1]
    if not df_stats.empty:
        df_plot = df_stats.copy()
        df_plot['-log10(q)'] = -np.log10(df_plot['q_value'].clip(lower=1e-10))
        
        colors = ['red' if q < 0.05 else 'gray' for q in df_plot['q_value']]
        
        ax2.scatter(df_plot['cliffs_delta'], df_plot['-log10(q)'], 
                   c=colors, alpha=0.6, s=50)
        ax2.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='q=0.05')
        ax2.axvline(x=0, color='black', linestyle='-', alpha=0.3)
        
        # Label significant points
        sig_points = df_plot[df_plot['q_value'] < 0.05]
        for _, row in sig_points.iterrows():
            if abs(row['cliffs_delta']) > 0.2:
                ax2.annotate(row['descriptor'], 
                           (row['cliffs_delta'], row['-log10(q)']),
                           fontsize=8, alpha=0.7)
        
        ax2.set_xlabel("Cliff's Delta\n← Associated with Poor Performance | Good Performance →")
        ax2.set_ylabel("-log10(q-value)")
        ax2.set_title(f'Descriptor Enrichment: {method_name}')
        ax2.legend()
    
    # 3. Top descriptors for good performance
    ax3 = axes[1, 0]
    if not df_stats.empty:
        # Get descriptors associated with good performance
        good_descriptors = df_stats[
            (df_stats['cliffs_delta'] > 0) & (df_stats['q_value'] < 0.1)
        ].sort_values('q_value').head(10)
        
        if not good_descriptors.empty:
            y_pos = np.arange(len(good_descriptors))
            ax3.barh(y_pos, good_descriptors['cliffs_delta'], color='green', alpha=0.7)
            ax3.set_yticks(y_pos)
            ax3.set_yticklabels(good_descriptors['descriptor'], fontsize=9)
            ax3.set_xlabel("Cliff's Delta (Effect Size)")
            ax3.set_title(f'Features Associated with Good {method_name} Performance')
        else:
            ax3.text(0.5, 0.5, 'No significant features for good performance',
                    ha='center', va='center', transform=ax3.transAxes)
    
    # 4. Top descriptors for poor performance
    ax4 = axes[1, 1]
    if not df_stats.empty:
        # Get descriptors associated with poor performance
        poor_descriptors = df_stats[
            (df_stats['cliffs_delta'] < 0) & (df_stats['q_value'] < 0.1)
        ].sort_values('q_value').head(10)
        
        if not poor_descriptors.empty:
            y_pos = np.arange(len(poor_descriptors))
            ax4.barh(y_pos, -poor_descriptors['cliffs_delta'], color='red', alpha=0.7)
            ax4.set_yticks(y_pos)
            ax4.set_yticklabels(poor_descriptors['descriptor'], fontsize=9)
            ax4.set_xlabel("Cliff's Delta (Effect Size, absolute)")
            ax4.set_title(f'Features Associated with Poor {method_name} Performance')
        else:
            ax4.text(0.5, 0.5, 'No significant features for poor performance',
                    ha='center', va='center', transform=ax4.transAxes)
    
    plt.suptitle(f'{method_name} Performance Analysis by Molecular Descriptors', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Use consistent naming based on method, not display name
    method_key = 'c3p' if 'C3P' in method_name else 'chebifier'
    plot_path = out_dir / f"{method_key}_performance_analysis.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    logger.info(f"Saved performance plots to {plot_path}")
    plt.close()


def generate_method_interpretation(df_cls: pd.DataFrame, df_stats: pd.DataFrame, 
                                  method_name: str, f1_threshold: float, out_dir: Path) -> None:
    """Generate interpretation for single method performance analysis."""
    summary_lines = []
    summary_lines.append("=" * 80)
    summary_lines.append(f"{method_name.upper()} PERFORMANCE ANALYSIS")
    summary_lines.append("=" * 80)
    summary_lines.append("")
    
    # Performance summary
    summary_lines.append("PERFORMANCE OVERVIEW:")
    summary_lines.append("-" * 40)
    
    perf_counts = df_cls["performance_group"].value_counts()
    total_classes = len(df_cls)
    
    summary_lines.append(f"Total classes analyzed: {total_classes}")
    summary_lines.append(f"Performance breakdown:")
    for group in ["good", "moderate", "poor"]:
        if group in perf_counts.index:
            count = perf_counts[group]
            pct = (count / total_classes) * 100
            summary_lines.append(f"  - {group.capitalize()} (F1 {'>' if group == 'good' else '<'} {f1_threshold if group != 'moderate' else '0.5-0.7'}): {count} classes ({pct:.1f}%)")
    
    summary_lines.append(f"\nMean F1 score: {df_cls['f1_score'].mean():.3f}")
    summary_lines.append(f"Median F1 score: {df_cls['f1_score'].median():.3f}")
    
    # Key molecular features
    summary_lines.append("\n\nKEY MOLECULAR FEATURES:")
    summary_lines.append("-" * 40)
    
    if not df_stats.empty:
        sig_good = df_stats[(df_stats['q_value'] < 0.05) & (df_stats['cliffs_delta'] > 0)].sort_values('q_value')
        sig_poor = df_stats[(df_stats['q_value'] < 0.05) & (df_stats['cliffs_delta'] < 0)].sort_values('q_value')
        
        if not sig_good.empty:
            summary_lines.append(f"\n✅ {method_name} PERFORMS WELL WITH:")
            for _, row in sig_good.head(5).iterrows():
                summary_lines.append(f"  • Higher {row['descriptor']}")
                summary_lines.append(f"    - Median in good classes: {row['median_good_perf']:.2f}")
                summary_lines.append(f"    - Median in poor classes: {row['median_poor_perf']:.2f}")
                summary_lines.append(f"    - Effect size: δ={row['cliffs_delta']:.3f}, q={row['q_value']:.4f}")
        
        if not sig_poor.empty:
            summary_lines.append(f"\n❌ {method_name} STRUGGLES WITH:")
            for _, row in sig_poor.head(5).iterrows():
                summary_lines.append(f"  • Higher {row['descriptor']}")
                summary_lines.append(f"    - Median in poor classes: {row['median_poor_perf']:.2f}")
                summary_lines.append(f"    - Median in good classes: {row['median_good_perf']:.2f}")
                summary_lines.append(f"    - Effect size: δ={row['cliffs_delta']:.3f}, q={row['q_value']:.4f}")
        
        if sig_good.empty and sig_poor.empty:
            summary_lines.append("\n⚠ No molecular features showed significant association with performance after FDR correction")
    
    # Write summary
    summary_text = "\n".join(summary_lines)
    # Use consistent naming based on method, not display name
    method_key = 'c3p' if 'C3P' in method_name else 'chebifier'
    summary_path = out_dir / f"{method_key}_interpretation.txt"
    with open(summary_path, 'w') as f:
        f.write(summary_text)
    
    logger.info(f"Saved interpretation to {summary_path}")
    
    # Print to console
    typer.echo(f"\n{method_name} Performance Summary:")
    typer.echo("-" * 40)
    typer.echo(f"Good performance (F1>{f1_threshold}): {perf_counts.get('good', 0)} classes")
    typer.echo(f"Poor performance (F1<0.5): {perf_counts.get('poor', 0)} classes")
    
    if not df_stats.empty:
        sig_features = df_stats[df_stats['q_value'] < 0.05]
        if not sig_features.empty:
            typer.echo(f"\nFound {len(sig_features)} significant molecular features")


@app.command()
def explore_differences(
    comparison_file: Path = typer.Argument(..., help="Path to pairwise_comparison.csv from c3p-compare"),
    dataset: Optional[Path] = typer.Option(None, "--dataset", "-d", help="Path to c3p dataset JSON file"),
    molecules: Optional[Path] = typer.Option(None, "--molecules", "-m", help="Alternative: Path to molecules CSV"),
    delta: float = typer.Option(0.05, help="|difference| threshold for favored groups"),
    agg: str = typer.Option("median", help="Class-level aggregation", case_sensitive=False),
    outdir: Path = typer.Option(Path("results/analysis"), help="Output directory"),
    include_metadata: bool = typer.Option(False, "--include-metadata", help="Include class metadata in analysis"),
    max_molecules_per_class: Optional[int] = typer.Option(None, "--max-molecules", help="Max molecules per class (for faster testing)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable debug logging"),
) -> None:
    """
    Analyze differences between ensemble and chebifier performance using c3p comparison output.
    
    This command simplifies the workflow by directly using pairwise_comparison.csv files
    and can automatically fetch ChEBI class data from c3p datasets.
    
    Examples:
        # Using dataset file
        python analyze_differences.py explore-differences \\
            results/2025/comparison-chebifier/pairwise_comparison.csv \\
            --dataset results/2024/inputs/dataset.json
        
        # Using molecules CSV
        python analyze_differences.py explore-differences \\
            results/2025/comparison-chebifier/pairwise_comparison.csv \\
            --molecules data/molecules.csv
    """
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="[%(levelname)s] %(message)s")
    
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Read the comparison file (pairwise_comparison.csv)
    df_perf = pd.read_csv(comparison_file)
    
    # Handle different possible column names
    if "chemical_class" not in df_perf.columns and "class" in df_perf.columns:
        df_perf = df_perf.rename(columns={"class": "chemical_class"})
    
    # Get molecules data
    df_mol = None
    if dataset and dataset.exists():
        logger.info(f"Loading dataset from {dataset}")
        ds = load_dataset(dataset)
        if ds:
            df_mol = dataset_to_molecules_df(ds)
            logger.info(f"Loaded {len(df_mol)} molecule-class pairs from dataset")
            
            # If include_metadata, add class metadata
            if include_metadata and ds:
                metadata_rows = []
                for cls in ds.get("classes", []):
                    metadata_rows.append({
                        "chemical_class": cls.get("name"),
                        "chebi_id": cls.get("id"),
                        "definition": cls.get("definition")
                    })
                df_metadata = pd.DataFrame(metadata_rows)
                df_perf = df_perf.merge(df_metadata, on="chemical_class", how="left")
    
    elif molecules and molecules.exists():
        logger.info(f"Loading molecules from {molecules}")
        df_mol = pd.read_csv(molecules)
    else:
        # Try to find dataset files automatically
        search_paths = [
            Path("results/2025/benchmark/dataset.json"),
        ]
        for path in search_paths:
            if path.exists():
                logger.info(f"Found dataset at {path}")
                ds = load_dataset(path)
                if ds:
                    df_mol = dataset_to_molecules_df(ds)
                    break
    
    if df_mol is None or df_mol.empty:
        logger.error("No molecule data found. Please provide --dataset or --molecules")
        raise typer.Exit(1)
    
    # Detect class column
    class_col = detect_class_column(df_mol)
    
    # Limit molecules per class if requested
    if max_molecules_per_class:
        logger.info(f"Limiting to {max_molecules_per_class} molecules per class")
        # Get classes from performance data
        perf_classes = set(df_perf["chemical_class"].unique())
        # Filter molecules to only include those classes
        df_mol = df_mol[df_mol[class_col].isin(perf_classes)]
        # Sample molecules per class
        sampled_rows = []
        for cls in perf_classes:
            cls_rows = df_mol[df_mol[class_col] == cls]
            if len(cls_rows) > max_molecules_per_class:
                cls_rows = cls_rows.sample(n=max_molecules_per_class, random_state=42)
            sampled_rows.append(cls_rows)
        if sampled_rows:
            df_mol = pd.concat(sampled_rows, ignore_index=True)
        logger.info(f"Reduced to {len(df_mol)} molecules total")
    
    # Build descriptors
    logger.info("Computing molecular descriptors...")
    df_desc = build_descriptor_table(df_mol)
    
    # Expand classes
    df_expanded = expand_classes(df_mol, class_col)
    
    # Aggregate at class level
    df_cls = aggregate_class_descriptors(df_expanded, df_desc, agg.lower())
    
    # Label groups based on performance differences
    df_labeled = label_groups(df_cls, df_perf, delta)
    
    # Save class-level descriptor table
    output_file = outdir / "class_descriptor_table.csv"
    df_labeled.to_csv(output_file, index=False)
    logger.info(f"Saved class descriptor table to {output_file}")
    
    # Run statistical tests
    df_stats = run_mwut(df_labeled, outdir)
    if not df_stats.empty:
        logger.info(f"Saved Mann-Whitney U test results to {outdir / 'mwut_enrichment.csv'}")
        
        # Create visualization plots
        if PLOTTING_AVAILABLE:
            create_enrichment_plots(df_stats, outdir)
        
        # Print top enriched descriptors
        typer.echo("\nTop enriched descriptors (q < 0.05):")
        sig_descriptors = df_stats[df_stats["q_value"] < 0.05]
        if not sig_descriptors.empty:
            for _, row in sig_descriptors.head(10).iterrows():
                direction = "↑ ensemble" if row["median_ensemble"] > row["median_chebifier"] else "↑ chebifier"
                typer.echo(f"  {row['descriptor']:25s} q={row['q_value']:.3f} {direction}")
        else:
            typer.echo("  No significantly enriched descriptors found")
    
    # Run logistic regression
    run_logreg(df_labeled, outdir)
    logger.info(f"Saved logistic regression coefficients to {outdir / 'logreg_coefficients.csv'}")
    
    # Run continuous regression analysis
    regression_results = run_continuous_regression(df_labeled, outdir)
    if not regression_results.empty:
        typer.echo("\nTop features from continuous regression:")
        for _, row in regression_results.head(10).iterrows():
            if row['abs_coefficient'] > 0.01:
                direction = "→ ensemble" if row['coefficient'] > 0 else "→ chebifier"
                typer.echo(f"  {row['descriptor']:30s} coef={row['coefficient']:+.3f} {direction}")
    
    # Generate interpretation summary
    generate_interpretation_summary(df_stats, regression_results, outdir)
    
    typer.echo(f"\n✓ Analysis complete. Results saved to: {outdir.resolve()}")


if __name__ == "__main__":
    app()
