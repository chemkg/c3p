import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import List, Optional, Annotated, Dict
from unittest.mock import right

import pandas as pd
import typer
from matplotlib import pyplot as plt
import seaborn as sns
from diskcache import Cache
from oaklib.datamodels.vocabulary import IS_A

from c3p.cli import configure_logging, verbose_option
from c3p.datamodel import Dataset, Config, ResultSet, EvaluationResult, CodeStatistics, Result, ChemicalStructure
from c3p.dumper import result_as_dict
from c3p.learn import safe_name, run_code
from c3p.outcome_stats import calculate_grouped_stats, assign_ranks_per_group
from c3p.plotting import create_scatter_matrix, plot_scatter, create_scatter_matrix_train_colored, plot_method_comparison_scatter

app = typer.Typer(help="Compare models.")

import logging

# Set up logger
logger = logging.getLogger()  # Get root logger

app = typer.Typer()


# Create a persistent cache for ChEBI metadata
chebi_cache = Cache('.chebi_metadata')

@chebi_cache.memoize()
def get_chebi_metadata(chebi_id: str) -> Dict[str, str]:
    """
    Get ChEBI metadata for a given ID, with automatic disk caching.
    
    Args:
        chebi_id: ChEBI identifier (e.g., 'CHEBI:12345')
        
    Returns:
        Dictionary with metadata fields
    """
    metadata_dict = {
        'definition': '',
        'formula': '',
        'xrefs': '',
        'curation_status': ''
    }

    print("Getting adapter..")
    try:
        from oaklib import get_adapter
        # Initialize adapter if not already done
        if not hasattr(get_chebi_metadata, '_adapter'):
            logger.info(f"Initializing ChEBI adapter (one-time setup, may take a moment)...")
            get_chebi_metadata._adapter = get_adapter("sqlite:obo:chebi")
            logger.info("ChEBI adapter initialized")
        
        adapter = get_chebi_metadata._adapter
        print(f"Looking up metadata for {chebi_id}")
        metadata = adapter.entity_metadata_map(chebi_id)
        
        metadata_dict['definition'] = metadata.get('IAO:0000115', [''])[0] if 'IAO:0000115' in metadata else ''
        metadata_dict['formula'] = metadata.get('obo:chebi/formula', [''])[0] if 'obo:chebi/formula' in metadata else ''
        metadata_dict['xrefs'] = '|'.join(metadata.get('oio:hasDbXref', []))
        metadata_dict['curation_status'] = '|'.join(metadata.get('oio:inSubset', []))
        print(f"Metadata for {chebi_id}: {metadata_dict}")
        
    except Exception as e:
        logger.debug(f"Could not get metadata for {chebi_id}: {e}")
    
    return metadata_dict


@app.command()
def combine(
        dirs: Optional[List[str]] = typer.Argument(..., help="directory name"),
        output_dir: Path = typer.Option(None, "--output", "-o", help="output directory"),
        exclude_incomplete: Optional[bool] = typer.Option(False, "--exclude-incomplete/--no-exclude-incomplete",
                                                          "-x/--no-x",
                                                          help="if set, exclude any result for a chemical not present in all sets"),
        verbose: Annotated[int, verbose_option] = 0
) -> None:
    """
    Evaluate a model on a dataset using a single class.
    """
    configure_logging(verbose)
    output_dir.mkdir(exist_ok=True, parents=True)
    rows = []
    for d in dirs:
        p = Path(d) / "eval"
        n = 0
        for fn in p.glob("*.json"):
            er = EvaluationResult.model_validate_json(fn.read_text())
            experiment_name = Path(d).name
            row = result_as_dict(er)
            row["experiment_name"] = experiment_name
            if er.train_results.best_result:
                row["chemical_name"] = er.train_results.best_result.chemical_class.name
            else:
                # chebifier results are test-only
                row["chemical_name"] = er.test_result.chemical_class.name
            if er.test_result.code:
                if not er.test_result.code_statistics:
                    er.test_result.code_statistics = CodeStatistics.from_code(er.test_result.code)
                for k, v in er.test_result.code_statistics.model_dump().items():
                    if isinstance(v, (int, float)):
                        row[k] = v
            rows.append(row)
            n += 1
        logger.info(f"Processed {n} results from {d}")
        if not n:
            logger.warning(f"No results found in {d}")
    if not rows:
        raise ValueError(f"No results found")


    df = pd.DataFrame(rows)
    df.to_csv(output_dir / "combined_results.csv", index=False)
    df.to_excel(output_dir / "combined_results.excel", index=False)

    pairwise_data = []
    for chem_class in df["chemical_class"].unique():
        subset = df[df["chemical_class"] == chem_class]
        for (left, right) in combinations(subset.itertuples(index=False), 2):
            pairwise_data.append({
                "chemical_class": chem_class,
                "left_method": left.experiment_name,
                "right_method": right.experiment_name,
                "f1_left": left.f1,
                "f1_right": right.f1,
                "difference": left.f1 - right.f1
            })

    # if a chemical count is < num_expts, exclude it
    expts = list(df["experiment_name"].unique())
    num_expts = len(expts)
    if exclude_incomplete:
        df = df.groupby("chemical_name").filter(lambda x: len(x) == num_expts)
        print(f"Filtering to {len(df)} rows with complete data")
    all_cls = list(df["chemical_class"].unique())

    # Create the pairwise comparison dataframe
    pairwise_df = pd.DataFrame(pairwise_data)
    pairwise_df.to_csv(output_dir / "pairwise_comparison.csv", index=False)

    top_n_threshold = 10
    top_differences1 = pairwise_df.groupby(["left_method", "right_method"]).apply(
        lambda x: x.nlargest(top_n_threshold, "difference")
    ).reset_index(drop=True)
    top_differences2 = pairwise_df.groupby(["left_method", "right_method"]).apply(
        lambda x: x.nsmallest(top_n_threshold, "difference")
    ).reset_index(drop=True)
    top_differences = pd.concat([top_differences1, top_differences2])
    top_differences.to_csv(output_dir / "top_pairwise_differences.csv", index=False)

    fig = create_scatter_matrix(df)
    fig.savefig(output_dir / "scatter_matrix.png")

    # Create train F1 colored scatter matrices for each experiment
    for experiment_name in df["experiment_name"].unique():
        # Skip experiments that don't have train_f1 data (like chebifier)
        exp_df = df[df["experiment_name"] == experiment_name]
        if "train_f1" not in exp_df.columns or exp_df["train_f1"].isna().all():
            logger.info(f"Skipping {experiment_name} - no train_f1 data available")
            continue
        
        fig = create_scatter_matrix_train_colored(df, experiment_name)
        if fig is not None:
            safe_exp_name = experiment_name.replace("/", "_").replace(" ", "_")
            fig.savefig(output_dir / f"scatter_matrix_{safe_exp_name}_train_scaled.png")
            plt.close(fig)
            logger.info(f"Created train F1 colored scatter matrix for {experiment_name}")
        else:
            logger.warning(f"Could not create scatter matrix for {experiment_name}")

    # Create pairwise method comparison scatter plots
    methods = list(df["experiment_name"].unique())
    for i, method1 in enumerate(methods):
        for j, method2 in enumerate(methods):
            if i < j:  # Only create plots for unique pairs
                fig = plot_method_comparison_scatter(df, method1, method2)
                if fig is not None:
                    safe_method1 = method1.replace("/", "_").replace(" ", "_")
                    safe_method2 = method2.replace("/", "_").replace(" ", "_")
                    fig.savefig(output_dir / f"method_comparison_{safe_method1}_vs_{safe_method2}.png")
                    plt.close(fig)
                    logger.info(f"Created method comparison plot: {method1} vs {method2}")
                else:
                    logger.warning(f"Could not create comparison plot for {method1} vs {method2}")

    # ranking

    for tt in ["test", "train"]:

        metric = "f1" if tt == "test" else "train_f1"
        rank_col = "rank" if tt == "test" else "train_rank"

        expt_stats = df.groupby('experiment_name').agg({
            metric: ['mean', 'std', 'min', 'max']
        }).round(2)
        expt_stats_by_mean = expt_stats.sort_values((metric, 'mean'))

        plt.figure(figsize=(8, 12))
        sns.violinplot(data=df, y=metric, x="experiment_name", order=expt_stats_by_mean.index, inner="box")
        plt.title(f'Distribution of {metric} by Experiment')
        plt.xticks(rotation=90, ha="center")  # `ha="center"` ensures text stays aligned
        plt.tight_layout()
        plt.savefig(output_dir / f"{metric}_distribution.png")

        assign_ranks_per_group(df, "experiment_name", metric=metric, col_name=rank_col)
        rank_stats = df.groupby('chemical_class').agg({
            rank_col: ['mean', 'std', 'min', 'max']
        }).round(2)

        df_sorted = df.sort_values('rank')
        df.to_csv(output_dir / f"combined_results_with_{rank_col}.csv", index=False)

        all_cls = list(df_sorted["chemical_class"].unique())
        num_cls = len(all_cls)


        # Add range column
        rank_stats[rank_col, 'range'] = rank_stats[rank_col, 'max'] - rank_stats[rank_col, 'min']

        rank_stats_by_std = rank_stats.sort_values((rank_col, 'std'))
        rank_stats_by_mean = rank_stats.sort_values((rank_col, 'mean'))
        print(f"Rank stats by std ({tt})")
        print(rank_stats_by_std)
        print(f"Rank stats by mean ({tt})")
        print(rank_stats_by_mean)
        with open(output_dir / f"{rank_col}_stats_by_std.csv", "w") as f:
            rank_stats_by_std.to_csv(f)
        with open(output_dir / f"{rank_col}_stats_by_mean.csv", "w") as f:
            rank_stats_by_mean.to_csv(f)

        plt.figure(figsize=(8, 50))
        sns.boxplot(data=df_sorted, y='chemical_class', x='rank', order=rank_stats_by_mean.index)
        #plt.xticks(rotation=90, ha='right')
        plt.title(f'Rank ({tt}) Distribution by Chemical Class, sorted by mean rank')
        plt.tight_layout()
        plt.savefig(output_dir / f"{rank_col}_distribution.png")

        top_n_variable = rank_stats.nsmallest(50, (rank_col, 'mean')).index
        df_subset = df[df['chemical_class'].isin(top_n_variable)]
        plt.figure(figsize=(8, 12))
        sns.boxplot(data=df_subset, y='chemical_class', x=rank_col, order=top_n_variable)
        # plt.xticks(rotation=90, ha='right')
        plt.title(f'Rank ({tt}) Distribution, top 50')
        plt.tight_layout()
        plt.savefig(output_dir / f"{rank_col}_distribution_top_50.png")

        plt.figure(figsize=(8, 50))
        sns.boxplot(data=df, y='chemical_class', x=rank_col, order=rank_stats_by_std.index)
        # plt.xticks(rotation=90, ha='right')
        plt.title(f'Rank ({tt}) Distribution by Chemical Class, sorted by std')
        plt.tight_layout()
        plt.savefig(output_dir / f"{rank_col}_distribution_by_std.png")

        # Print summary stats
        print(f"Chemical classes with most consistent rankings (lowest std) ({tt}):")
        print(rank_stats.nsmallest(5, (rank_col, 'std')))

        print(f"\nChemical classes with most variable rankings (highest std) ({tt}):")
        print(rank_stats.nlargest(5, (rank_col, 'std')))

        # Optional: Calculate correlation between mean rank and rank variability
        correlation = rank_stats[rank_col]['mean'].corr(rank_stats[rank_col]['std'])
        print(f"\nCorrelation between mean {rank_col} and rank std: {correlation:.3f}")

    df_sorted = df.sort_values('complexity')
    plt.clf()
    plt.figure(figsize=(8, 50))
    sns.boxplot(data=df_sorted, y='chemical_class', x='complexity')
    plt.title('Code Complexity Distribution by Class')
    plt.tight_layout()
    plt.savefig(output_dir / "complexity_distribution.png")

    #rank_stats = rank_stats.reset_index()
    #rank_stats = rank_stats.sort_values(('std',))
    #print("Rank stats sorted by rank std")
    #print(rank_stats)
    print(df["experiment_name"].unique())
    stats = calculate_grouped_stats(df, "experiment_name")
    print(stats)
    print(stats.T)
    stats.to_csv(output_dir / "summary_stats.csv")
    stats.to_excel(output_dir / "summary_stats.xlsx")
    cols = ['micro_f1', 'macro_f1', 'micro_precision', 'macro_precision', 'micro_recall', 'macro_recall', 'micro_accuracy', 'macro_accuracy', 'micro_negative_predictive_value', 'macro_negative_predictive_value']
    stats_slim = stats[cols]
    stats_slim.to_csv(output_dir / "summary_stats_slim.csv")
    stats_slim.to_excel(output_dir / "summary_stats_slim.xlsx")

    stats.reset_index(inplace=True)
    stats.rename(columns={"index": "model"}, inplace=True)
    # rename values in model column, e.g. gpt-4o => 4o

    def rename_model(this_df):
        this_df["model"] = this_df["model"].str.replace("gpt-4o", "4o").str.replace("-preview", "").str.replace("-undef", "").str.replace("ensembl-5", "ensemble").str.replace("-hi", "-iter6").str.replace("-coder", "").str.replace("-3-sonnet", "")

    #stats["model"] = stats["model"].str.replace("gpt-4o", "4o").str.replace("-preview", "").str.replace("-undef", "").str.replace("ensembl-5", "ensemble").str.replace("-hi", "-iter6").str.replace("-coder", "").str.replace("-3-sonnet", "")
    rename_model(stats)
    #all_cls = list(stats["chemical_class"].unique())
    #num_cls = len(all_cls)
    for typ in ["micro", "macro"]:
        main_metric = f"{typ}_f1"
        # Sort dataframe by 'f1' in descending order
        df_sorted = stats.sort_values(by=main_metric, ascending=False).reset_index(drop=True)

        # Plotting sorted bar chart
        fig, ax = plt.subplots(figsize=(10, 6))
        width = 0.25  # Bar width
        x = range(len(df_sorted))

        # Define distinct colors for the metrics
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

        # Plot each metric as a separate group of bars with distinct colors
        ax.bar([i - width for i in x], df_sorted[main_metric], width, label=f'{typ} f1', color=colors[0])
        ax.bar(x, df_sorted[f'{typ}_precision'], width, label=f'{typ} precision', color=colors[1])
        ax.bar([i + width for i in x], df_sorted[f'{typ}_recall'], width, label=f'{typ} recall', color=colors[2])

        # Add labels and title
        ax.set_xlabel("Models")
        ax.set_ylabel("Scores")

        ax.set_title(f"Comparison using {typ} stats, sorted by F1")
        ax.set_xticks(x)
        ax.set_xticklabels(df_sorted['model'])
        ax.legend()

        plt.xticks(rotation=90, ha="center")  # `ha="center"` ensures text stays aligned

        # Display the plot
        plt.tight_layout()
        #plt.show()
        plt.savefig(output_dir / f"model_comparison_{typ}.png")

    # Code Complexity plot
    # Melt the DataFrame to have a "long" format for sns.barplot
    df["model"] = df["experiment_name"]
    rename_model(df)
    complexity_metrics = ['complexity', "max_indent",  'log_lines_of_code', 'methods_called_count', 'returns_count', 'smarts_strings_count']
    model_order = df.groupby('model')['complexity'].mean().sort_values(ascending=False).index
    df_melted = df.melt(id_vars='model', value_vars=complexity_metrics,
                        var_name='metric', value_name='value')

    # Plotting
    plt.figure(figsize=(10, 6))
    sns.barplot(x='model', y='value', hue='metric', data=df_melted, order=model_order, palette='viridis')
    plt.title('Code Metrics per Model')
    plt.xlabel('Model')
    plt.ylabel('Value')
    plt.legend(title='Metric')
    plt.xticks(rotation=90, ha='center')
    #plt.show()
    plt.savefig(output_dir / f"complexity_comparison.png")

    # save 50 most complex
    most_complex = df.sort_values('complexity', ascending=False).head(50)[["model", "chemical_class"] + complexity_metrics]
    most_complex.to_csv(output_dir / "most_complex.csv", index=False)

    # least complex
    least_complex = df.sort_values('complexity', ascending=True).head(50)[["model", "chemical_class"] + complexity_metrics]
    least_complex.to_csv(output_dir / "least_complex.csv", index=False)

    for m in complexity_metrics:
        plt.figure(figsize=(6, 6))
        plot_scatter(df, m, "f1", results_dir=output_dir)


if __name__ == "__main__":
    app()