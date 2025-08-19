import os
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Annotated, Dict
import json

import pandas as pd
import typer
import numpy as np

from c3p.cli import verbose_option, configure_logging
from c3p.datamodel import Dataset, Config, EvaluationResult, ResultSet, Result
from c3p.dumper import result_as_dict, write_programs, write_programs_for_commits
from c3p.learn import learn_ontology_iter, evaluate_class

app = typer.Typer(help="CHEBI learn.")

import logging

# Set up logger
logger = logging.getLogger()  # Get root logger


app = typer.Typer()


@app.command()
def validate(
        working_dir: Path = typer.Option(None, "--workdir", "-w", help="path to workdir"),
        dataset_path: Path = typer.Option(None, "--dataset", "-d", help="path to dataset"),
        verbose: Annotated[int, verbose_option] = 0
) -> None:
    """
    Evaluate a model on a dataset using a single class.
    """
    n = 0
    eval_dir = working_dir.parent / "eval"
    eval_dir.mkdir(exist_ok=True, parents=True)
    configure_logging(verbose)
    eval_results = []
    if dataset_path:
        with open(dataset_path, "r") as f:
            dataset = Dataset.model_validate_json(f.read())
            print(f"Classes: {len(dataset.classes)} Instances: {len(dataset.structures)}")
    ers = []
    for fn in working_dir.glob("*.json"):
        n = Path(fn).name
        ofn = eval_dir / n
        if ofn.exists():
            logger.info(f"Reusing {ofn}")
            er = EvaluationResult.model_validate_json(ofn.read_text())
        else:
            rset = ResultSet.model_validate_json(fn.read_text())
            er = evaluate_class(rset, dataset)
            ofn.write_text(er.model_dump_json(indent=2))
        ers.append(er)
        er_dict = result_as_dict(er)
        timestamp = os.path.getmtime(fn)
        # Convert to datetime object
        dt = datetime.fromtimestamp(timestamp)
        formatted_time = dt.strftime("%Y-%m-%d %H:%M:%S")
        er_dict["time"] = formatted_time
        eval_results.append(er_dict)
    if not eval_results:
        raise ValueError("No results found")
    df = pd.DataFrame(eval_results)
    print(df.describe())
    write_programs(ers, working_dir.parent)
    write_programs_for_commits(ers, working_dir.parent)




if __name__ == "__main__":
    app()