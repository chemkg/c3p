import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Annotated

import pandas as pd
import typer
import yaml
from matplotlib import pyplot as plt
from oaklib import get_adapter

from c3p.cli import verbose_option, configure_logging
from c3p.datamodel import Dataset, Config, EvaluationResult, ResultSet, CodeStatistics
from c3p.dumper import result_as_dict
from c3p.learn import safe_name
from c3p.outcome_stats import calculate_all_outcome_stats
from c3p.plotting import plot_scatter

app = typer.Typer(help="CHEBI learn.")

import logging

# Set up logger
logger = logging.getLogger()  # Get root logger


app = typer.Typer()


@app.command()
def summarize(
        working_dir: Path = typer.Option(None, "--workdir", "-w", help="path to workdir"),
        output_dir: Optional[Path] = typer.Option(None, "--output", "-o", help="output directory"),
        verbose: Annotated[int, verbose_option] = 0
) -> None:
    """
    Evaluate a model on a dataset using a single class.
    """
    n = 0
    dataset = Dataset(**json.load(open(working_dir.parent.parent / "benchmark" / "dataset.json")))
    configure_logging(verbose)
    eval_results = []
    eval_objects = []
    for fn in working_dir.glob("*.json"):
        if "/cache" in str(working_dir):
            er = ResultSet.model_validate_json(fn.read_text())
        else:
            er = EvaluationResult.model_validate_json(fn.read_text())
        eval_objects.append(er)
        #auprc = calculate_auprc(er)
        er_dict = result_as_dict(er)
        if isinstance(er, EvaluationResult):
            br = er.train_results.best_result
            if not br and "chebifier" in str(working_dir):
                br = er.test_result
        else:
            br = er
        if not br.code_statistics and br.code:
            print(f"calculated code statistics for {fn}")
            br.code_statistics = CodeStatistics.from_code(br.code)
            for k, v in br.code_statistics.model_dump().items():
                if isinstance(v, (int, float)):
                    er_dict[k] = v
        er_dict["num_train_examples"] = br.num_true_positives + br.num_false_negatives
        #if not er_dict["num_train_examples"]:
        #    print(f"no train examples for {er_dict['f1']} {fn}")
        timestamp = os.path.getmtime(fn)
        # Convert to datetime object
        dt = datetime.fromtimestamp(timestamp)
        formatted_time = dt.strftime("%Y-%m-%d %H:%M:%S")
        er_dict["time"] = formatted_time
        #er_dict["auprc"] = auprc
        eval_results.append(er_dict)
    df = pd.DataFrame(eval_results)
    df["log_num_train_examples"] = df["num_train_examples"].apply(lambda x: 0 if x == 0 else math.log(x) / math.log(10))
    print(df.describe())
    metrics = ["f1", "precision", "recall", "accuracy", "negative_predictive_value", "train_f1"]
    print(df.aggregate({k: ["mean", "std", "max", "min"] for k in metrics}))
    all_stats_sf = calculate_all_outcome_stats(df)
    print(all_stats_sf.T.round(3))
    if output_dir:
        output_dir.mkdir(exist_ok=True, parents=True)
        md_dir = output_dir / "markdown"
        md_dir.mkdir(exist_ok=True, parents=True)
        for eo in eval_objects:
            br = eo.train_results.best_result
            md_fn = md_dir / f"{safe_name(br.chemical_class.name)}.md"
            with open(md_fn, "w") as f:
                f.write(eo.markdown)
        df.to_csv(output_dir / "results.csv", index=False)
        df.to_excel(output_dir / "results.xslx", index=False)
        all_stats_sf.to_csv(output_dir / "outcome_stats.csv", index=False)
        df.describe().to_csv(output_dir / "summary.csv", index=True)
        df.describe().to_html(output_dir / "summary.html", index=True)
        plot_scatter(df, "train_f1", "f1",
                     title="F1 in training vs F1 on test data",
                     xlabel="Train F1", ylabel="Test F1",
                     results_dir = output_dir,
                     )
        plot_scatter(df, "log_num_train_examples", "f1", results_dir=output_dir)
        plot_scatter(df, "complexity", "f1", results_dir=output_dir)
        #plot_scatter(df, "num_train_examples", "f1", results_dir=output_dir)
        for min_f1 in [0.05, 0.5]:
            filtered_df = df[df["train_f1"] > min_f1]
            filtered_df.to_csv(output_dir / f"results_{min_f1}.csv", index=False)
            filtered_df.to_excel(output_dir / f"results_{min_f1}.xslx", index=False)
            filtered_df.to_html(output_dir / f"results_{min_f1}.html", index=False)
            filtered_df.describe().to_csv(output_dir / f"summary_{min_f1}.csv", index=False)
            filtered_df.describe().to_html(output_dir / f"summary_{min_f1}.html", index=False)
        df_sorted = df.sort_values("f1", ascending=False)
        df_sorted.to_csv(output_dir / "results_sorted.csv", index=False)
        # show top N chemicals
        N = 100
        df_sorted[:N][["chemical_class", "f1"]].to_csv(output_dir / "top_chemicals.csv", index=False)
        df_sorted[:N][["chemical_class", "f1", "false_positives", "false_negatives"]].to_csv(output_dir / "top_chemicals_failed.csv", index=False)
        df_sorted[:N][["chemical_class", "f1", "false_positives", "sample_false_negatives"]].to_html(output_dir / "top_chemicals_failed.html", index=False)
        # show bottom N chemicals
        df_sorted[-N:][["chemical_class", "f1"]].to_csv(output_dir / "bottom_chemicals.csv", index=False)
        df_sorted[["chemical_class", "f1"]].to_csv(output_dir / "ranked_chemicals.csv", index=False)

        df_sorted = df.sort_values("train_f1", ascending=False)
        # Define N (e.g., top and bottom 20)
        N = 30

        # Select the bottom N and top N rows
        df_bottom = df_sorted.head(N)
        df_top = df_sorted.tail(N)

        # Optional: add a column to distinguish between the two groups
        df_bottom = df_bottom.copy()  # Avoid SettingWithCopyWarning
        df_top = df_top.copy()
        df_bottom['Group'] = 'Bottom'
        df_top['Group'] = 'Top'

        # Combine the two DataFrames
        df_plot = pd.concat([df_bottom, df_top])

        # Plot a horizontal bar chart using Seaborn
        plt.figure(figsize=(10, 8))
        import seaborn as sns
        sns.barplot(
            data=df_plot,
            x='train_f1',
            y='chemical_class',
            hue='Group',  # Color-code to distinguish top vs. bottom
            dodge=False  # Keep the bars aligned on the same axis
        )
        plt.xlabel('Train F1 Score')
        plt.ylabel('Chemical Class')
        plt.title(f'Top and Bottom {N} Chemical Classes by Train F1')
        plt.tight_layout()
        #plt.show()
        plt.savefig(output_dir / "train_f1.png")

        rows = []
        for _, row in df_sorted.iterrows():
            rows.append(add_image_urls(row, dataset))

        img_df = pd.DataFrame(rows)

        N = 10
        for subset in ["all", "top", "bottom"]:
            if subset == "all":
                subset_df = img_df
            elif subset == "top":
                subset_df = img_df.head(N)
            else:
                # bottom, but also sort such that the worst is at the top
                subset_df = img_df.tail(N).sort_values("f1")

            # Apply formatting
            styled_df = subset_df.style.format({})

            if subset == "top":
                # green background for top rows
                color = "#4CAF50"
            elif subset == "bottom":
                # red background for bottom rows
                color = "#f44336"
            else:
                color = "#f0f0f0"

            styled_df.set_table_styles([
                {"selector": "th", "props": [("background-color", f"{color}"), ("color", "white"), ("padding", "10px")]},
                {"selector": "td", "props": [("word-wrap", "break-word"), ("white-space", "normal"), ("max-width", "250px"), ("font-size", "14px"), ("vertical-align", "middle")]},  # Align content
                #{"selector": "td", "props": [("padding", "10px"), ("border", "1px solid #ddd"), ("text-align", "center")]},
                {"selector": "table", "props": [("border-collapse", "collapse"), ("width", "100%")]},
            ])

            # Save to HTML
            html_output = styled_df.to_html(output_dir / f"classes_{subset}.html", escape=False, render_links=True, index=False)


chebi = get_adapter("sqlite:obo:chebi")

IMG_URL_TMPL = "https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId={n}"

def img_url(chebi_id: str) -> str:
    url = IMG_URL_TMPL.format(n=chebi_id.replace("CHEBI:", ""))
    return f'<img src="{url}" width="120" style="max-height: 100px; object-fit: contain; display: block; margin: auto;" onerror="this.style.display=\'none\';"/>'


def cell_div(text: str, img: str) -> str:
    return f'''
        <div style="display: flex; flex-direction: column; align-items: center; text-align: center;">
            <span style="font-size: 12px; margin-bottom: 5px;">{text}</span>
            {img}
        </div>
        '''

def add_image_urls(row: dict, dataset: Dataset):
    cc = dataset.get_chemical_class_by_name(row["chemical_class"])
    cc_id = cc.id
    example_smiles = cc.all_positive_examples[0]
    inst = dataset.smiles_to_instance()[example_smiles]
    chebi_ids = chebi.curies_by_label(inst.name)
    f1 = row["train_f1"]
    inst_img_href = img_url(chebi_ids[0])
    inst_name = inst.name
    # we will use precision 3 digits for f1

    new_row = {
        "f1": f"{f1:.3f}",
        "class": f"{cc.name}<br/>{img_url(cc_id)}",
        "example": f"{inst_name}<br/>{inst_img_href}",
    }
    new_row = {
        "f1": f"{f1:.3f}",
        "class": cell_div(cc.name, img_url(cc_id)),
        "example": cell_div(inst_name, inst_img_href),
    }
    return new_row




if __name__ == "__main__":
    app()