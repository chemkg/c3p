import json
import logging
import random
from copy import copy, deepcopy
from pathlib import Path
from typing import Tuple, Optional, Set, Dict, Any
from typing import List
import sys
import io
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from typing import Iterator
from pathlib import Path

from timeout_decorator import timeout

from c3p.curated import CURATED_PATH

test_dir = Path("tmp")

logger = logging.getLogger(__name__)

from c3p.datamodel import ChemicalStructure, Result, ChemicalClass, Config, ResultSet, EvaluationResult, Dataset, \
    SMILES_STRING, Outcome

example_dir = CURATED_PATH
VALIDATED_PROGRAM_EXAMPLES = ["triglyceride"]


BASE_SYSTEM_PROMPT = """
Write a program to classify chemical entities of a given class based on their SMILES string.

The program should consist of import statements (from rdkit) as well as a single function, named
is_<<chemical_class>> that takes a SMILES string as input and returns a boolean value plus a reason for the classification.

If the task is too hard or cannot be done, the program MAY return (None, None).

You are encouraged to think step by step. You can explain your reasoning BEFORE the start of the code
block. The code block should be enclosed by triple backticks (```) and MUST be in Python.
It must be a SINGLE valid code block. Do not include multiple code blocks, or split a code block.

You are also encouraged to include comments within the code.
"""

def safe_name(name: str) -> str:
    """
    Convert a name to a safe format for use in python functions.

    Example:

        >>> safe_name("foo' 3<->x bar")
        'foo__3___x_bar'

    """
    return "".join([c if c.isalnum() else "_" for c in name])

def generate_system_prompt(program_examples: List[str] = None):
    """
    Generate a system prompt for classifying chemical entities based on SMILES strings.

    Args:
        program_examples:

    Returns:
    """
    if program_examples is None:
        program_examples = VALIDATED_PROGRAM_EXAMPLES
    for example in program_examples:
        system_prompt = BASE_SYSTEM_PROMPT
        system_prompt += f"Here is an example for the chemical class {example}:\n{example}.py\n---\n"
        with open(f"{example_dir}/{example}.py", "r") as f:
            system_prompt += f.read()
    return system_prompt


def generate_main_prompt(chemical_class: str, definition: str, positive_instances: List[ChemicalStructure], err=None, prog=None):
    # replace all non-alphanumeric characters with underscores
    chemical_class_safe = safe_name(chemical_class)
    prompt = f"""
    Now create a program that classifies chemical entities of the class {chemical_class},
    defined as '{definition}'. The name of the function should be `is_{chemical_class_safe}`.
    Examples of structures that belong to this class are:
    """
    if definition is None:
        prompt = prompt.replace(", defined as 'None'", "")
    for instance in positive_instances:
        prompt += f" - {instance.name}: SMILES: {instance.smiles}\n"
    if err:
        prompt += f"\nYour last attempt failed with the following error: {err}\n"
    if prog:
        prompt += f"\nFor reference, this is the previous code:\n```\n{prog}\n```\n"
    return prompt

@timeout(2)
def eval_with_timeout(*args) -> Any:
    return eval(*args)

def run_code(code_str: str, function_name: str, args: List[Any], neg_args: List[Any], max_false: Optional[int]=None) -> List[Tuple[str, bool, str, dict]]:
    """
    Run the generated code and return the results.

    This expects the code to define a function with the name `function_name` that takes a single argument
    (SMILES) and returns a tuple of (arg, boolean, reason).

    Example:

        >>> code = "def is_foo(x): return x == 'foo', 'it is foo'"
        >>> run_code(code, 'is_foo', ['foo', 'bar'], [])
        [('foo', True, 'it is foo', {}), ('bar', False, 'it is foo', {})]

    Args:
        code_str:
        function_name:
        args:
        neg_args:
        max_false:

    Returns:

    """
    # suppress rdkit logging
    #log_stream = io.StringIO()
    #handler = logging.StreamHandler(log_stream)
    #tmp_logger = logging.getLogger()
    #tmp_logger.addHandler(handler)
    #tmp_logger.setLevel(logging.FATAL)
    # suppress at C++ level
    # https://github.com/rdkit/rdkit/issues/2683
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')

    __metadata__ = {}
    # redirect stdout to a StringIO object
    f = io.StringIO()
    with redirect_stdout(f), redirect_stderr(f):
        exec(code_str, globals())
    logger.info(f"Compilation Output: {f.getvalue()}")

    f = io.StringIO()
    with redirect_stdout(f), redirect_stderr(f):
        metadata = globals().get('__metadata__', {})
        vals = []
        num_false_positives = 0
        num_false_negatives = 0
        logger.info(f"Running code with {len(args)} + {len(neg_args)} instances")
        logger.debug(f"CODE:\n{code_str}")
        for i, arg in enumerate(args + neg_args):
            if max_false and i < len(args) and num_false_negatives > max_false:
                continue
            # use json.dumps to make sure that the argument is safely encoded
            # e.g. SMILES may use backslashes
            #arg_encoded = json.dumps(arg)
            arg_encoded = repr(arg)
            func_exec_str = f"{function_name}({arg_encoded})"
            logger.debug(f"Executing: {func_exec_str}")
            try:
                is_cls, reason = eval_with_timeout(func_exec_str)
            except Exception as e:
                raise RuntimeError(f"Error executing {func_exec_str}:\n {e}")
            if is_cls and i >= len(args):
                # in 2nd batch, should be negative
                num_false_positives += 1
                if max_false and num_false_positives > max_false:
                    break
            if not is_cls and i < len(args):
                # in 1st batch, should be positive
                num_false_negatives += 1

            vals.append((arg, is_cls, reason, metadata))
    logger.info(f"Exec Output: {f.getvalue()}")
    #tmp_logger.removeHandler(handler)
    return vals


@contextmanager
def capture_output():
    """
    Capture stdout and stderr using a context manager.

    Example:

        >>> with capture_output() as (out, err):
        ...     print("Hello, World!")
        ...
        >>> print(out.getvalue())
        Hello, World!

    """
    # Create StringIO objects to capture output
    stdout, stderr = io.StringIO(), io.StringIO()

    # Save the current stdout/stderr
    old_stdout, old_stderr = sys.stdout, sys.stderr

    try:
        # Replace stdout/stderr with our StringIO objects
        sys.stdout, sys.stderr = stdout, stderr
        yield stdout, stderr
    finally:
        # Restore the original stdout/stderr
        sys.stdout, sys.stderr = old_stdout, old_stderr

def evaluate_program(code_str: str, chemical_class: ChemicalClass, positive_instances: List[ChemicalStructure], negative_instances: List[ChemicalStructure], attempt=0,
                     threshold: Optional[float] = None,
                     model: Optional[Result] = None,
                     ) -> Result:
    """
    Evaluate a program on a set of positive and negative instances.

    Args:
        code_str:
        chemical_class:
        positive_instances:
        negative_instances:
        attempt: attempt number
        threshold: threshold for probability score, if None, no threshold is applied

    Returns:
        a Result object

    """
    if threshold == 0.0:
        threshold = None
    if threshold is None and not model:
        raise ValueError("Threshold is None, but no model is provided to get the threshold from")
    safe = safe_name(chemical_class.name)
    func_name = f"is_{safe}"
    positive_structures = [instance.smiles for instance in positive_instances]
    negative_structures = [instance.smiles for instance in negative_instances]
    smiles_to_instance = {instance.smiles: instance for instance in positive_instances + negative_instances}
    try:
        logger.info(f"Running code len: {len(code_str)} on {len(positive_structures)} + {len(negative_structures)} instances")
        with capture_output() as (stdout, stderr):
            #inputs = [instance.smiles for instance in positive_instances + negative_instances]
            results = run_code(code_str, func_name, positive_structures, negative_structures)
            logger.info(f"Results: {len(results)}")
    except Exception as e:
        # problem executing; we still yield a result, as this
        # may be useful for post-processing all results;
        # we also try again with a new attempt, unless we are suppressing LLM
        return Result(
            chemical_class=chemical_class,
            # config=config,
            code=code_str,
            # message=err,
            attempt=attempt,
            success=False,
            error=str(e) + stderr.getvalue(),
            stdout=stdout.getvalue(),
        )
    all_structures = positive_structures + negative_structures
    all_structures_count = len(all_structures)
    def mk(smiles, reason):
        instance = smiles_to_instance.get(smiles)
        return Outcome(smiles=smiles, name=instance.name, reason=reason)
    true_positives = [mk(smiles, reason) for smiles, is_cls, reason, _ in results if is_cls and smiles in positive_structures]
    true_negatives = [mk(smiles, reason) for smiles, is_cls, reason, _ in results if
                      not is_cls and not smiles in positive_structures]
    false_positives = [mk(smiles, reason) for smiles, is_cls, reason, _ in results if is_cls and not smiles in positive_structures]
    false_negatives = [mk(smiles, reason) for smiles, is_cls, reason, _ in results if not is_cls and smiles in positive_structures]


    # We avoid placing all negatives in the payload as these can be large
    result = Result(
        chemical_class=chemical_class,
        # config=config,
        code=code_str,
        # message=err,
        true_positives=true_positives,
        false_positives=false_positives,
        sample_false_negatives=false_negatives[:10],
        sample_true_negatives=true_negatives[:10],
        num_true_negatives=len(true_negatives),
        num_false_negatives=len(false_negatives),
        attempt=attempt,
        #stdout=stdout.getvalue(),
        error=stderr.getvalue(),
        success=True,
    )
    result.calculate()
    return result

def learn_program_single_iter(
        cls: ChemicalClass,
        positive_instances: List[ChemicalStructure],
        negative_instances: List[ChemicalStructure],
        attempt=0,
        err=None,
        prog=None,
        config: Optional[Config]=None,
) -> Iterator[Result]:
    """
    Main workflow

    The main workflow is a cycle between

    1. Generating code
    2. Running the code on positive and negative examples
    3. Go to 1 until `N` iterations or sufficient accuracy is received

    Each cycle will result a result. The final one is typically the best, but this
    is not guaranteed.

    :param cls: target chemical class for which is write a program
    :param positive_instances: positive instances
    :param negative_instances: negative instances
    :param attempt: counts which attempt this is
    :param err: error from previous iteration
    :param prog: program from previous iteration
    :param config: setup
    :return: iterator of results
    """
    logger.info(f"Test: {cls.name} attempt={attempt} err={err} prog={prog}")
    from llm import get_key, get_model
    if config is None:
        config = Config()
    cls_lite = cls.lite_copy()
    next_attempt = attempt + 1
    if next_attempt > config.max_attempts:
        print(f"FAILED: {cls.name} err={err[0:40]}")
        return

    if not positive_instances or not negative_instances:
        raise AssertionError("Run inference to get negative examples")
    system_prompt = generate_system_prompt(VALIDATED_PROGRAM_EXAMPLES)
    defn = cls.definition if config.use_definitions else None
    main_prompt = generate_main_prompt(cls.name, defn, positive_instances[0:config.max_positive_in_prompt], err=err, prog=prog)
    logger.info(f"System Prompt: {system_prompt}")
    logger.info(f"Main Prompt: {main_prompt}")
    model = get_model(config.llm_model_name)
    if model.needs_key:
        model.key = get_key(None, model.needs_key, model.key_env_var)
    options = {}
    if "llama" in config.llm_model_name:
        options = {"max_tokens": 3000}
    if "deepseek" in config.llm_model_name.lower():
        options = {"max_tokens": 8192 * 2}
    if "o1" in config.llm_model_name:
        response = model.prompt(f"SYSTEM PROMPT: {system_prompt} MAIN PROMPT: {main_prompt}", stream=False)
    else:
        response = model.prompt(main_prompt, system=system_prompt, **options)
    code_str = response.text()
    if not code_str:
        print(f"No code returned for {cls.name} // {response}")
    reasoning = None
    if "```" in code_str:  # Remove code block markdown
        parts = code_str.split("```")
        reasoning = parts[0].strip()
        code_str = parts[1].strip()
        if code_str.startswith("python"):
            code_str = code_str[6:]
        code_str = code_str.strip()

    code_str = repair_code(code_str)
    result = evaluate_program(code_str, cls_lite, positive_instances, negative_instances)
    result.attempt = attempt
    result.reasoning = reasoning
    if err:
        result.message = err
    logger.info(f"Attempt {attempt} for {cls.name} F1={result.f1}")
    yield result
    if result.f1 is None or result.f1 < config.f1_threshold:
        max_examples = config.max_examples_in_feedback
        # try again, feeding in the results of the current attempt
        msg = ""
        if result.error:
            msg += f"\nError: {result.error}"
        msg += f"\nAttempt failed: F1 score of {result.f1} is too low."
        msg += f"\nOutcomes:\n------\n"
        def ser_list(outcomes: List[Outcome], categ: str):
            if not outcomes:
                return "NONE"
            return "\n * ".join(f"SMILES: {o.smiles} NAME: {o.name} REASON: {categ} {o.reason}" for o in outcomes[:max_examples])

        msg += "\nTrue positives: " + ser_list(result.true_positives, "CORRECT")
        msg += "\nFalse positives: " + ser_list(result.false_positives, "WRONGLY CLASSIFIED")
        msg += "\nFalse negatives: " + ser_list(result.sample_false_negatives, "MISSED")
        msg += f"\n------\n"
        msg += "\nIn your reasoning step, analyze the previous program and the above outcomes, "
        msg += "hypothesizing about what went wrong, and how to improve.\n"
        if config.use_the_force:
            msg += "IMPORTANT NOTE: I do not have 100% confidence in the benchmark I am using. "
            msg += "There may be occasional and systematic mistakes. "
            msg += "Use your best judgment, and if you think the classifications your program are "
            msg += "consistent with your understanding if the meaning of the chemical class, then "
            msg += "you can ignore outliers, but explain your reasoning in doing so. "
            msg += "I have great confidence in your broad understanding of chemistry and your ability to "
            msg += "translate this into code."
        logger.info(f"Retrying {cls.name} with new prompt")
        yield from learn_program_single_iter(cls, positive_instances, negative_instances, config=config, attempt=next_attempt, err=msg, prog=code_str)


def learn_program(
        cls: ChemicalClass,
        positive_instances: List[ChemicalStructure],
        negative_instances: List[ChemicalStructure],
        config: Config) -> ResultSet:
    """
    Learn a program for a chemical class.

    Args:
        cls:
        positive_instances:
        negative_instances:
        config:

    Returns:

    """
    # print(f"## {test_cls.name} POS={len(test_cls.instances)} NEG={len(test_cls.negatives)}")
    results = []
    for result in learn_program_single_iter(cls, positive_instances, negative_instances, config=config):
        logger.info(f"attempt={result.attempt} compiled={result.success} tp={result.num_true_positives} tn={result.num_true_negatives} fp={result.num_false_positives} f1={result.f1}, len=={len(result.code)}")
        results.append(result)
        result.calculate()
    return ResultSet.from_results(results)

HALLUCINATED_CODE_LINES = [
    ("from rdkit.Chem import rdDecomposition", "rdCombination is hallucinated, it does not exist"),
]
def repair_code(code: str) -> str:
    """
    Repair code by adding comments for hallucinated code lines.

    Args:
        code:

    Returns:

    """
    for hcl, reason in HALLUCINATED_CODE_LINES:
        if hcl not in code:
            code = code.replace(hcl, f"# {hcl} ## {reason}")
    return code

def randomize_and_split_list(l: List, proportion: float) -> Tuple[List, List]:
    l = copy(l)
    random.shuffle(l)
    n = int(len(l) * proportion)
    return l[:n], l[n:]



def get_positive_and_negative_train_instances(cls: ChemicalClass, dataset: Dataset) -> Tuple[List[ChemicalStructure], List[ChemicalStructure]]:
    """
    Get positive and negative instances for a chemical class.

    Args:
        cls:
        dataset:

    Returns:

    """
    s2i = dataset.smiles_to_instance()
    all_validation = set(dataset.validation_examples)
    all_positive = set(cls.all_positive_examples)
    all_smiles = dataset.all_smiles()
    positive_examples = list(all_positive - all_validation)
    positive_instances = [s2i[smiles] for smiles in positive_examples]
    negative_examples = list((all_smiles - all_positive) - all_validation)
    negative_instances = [s2i[smiles] for smiles in negative_examples]
    return positive_instances, negative_instances


def get_positive_and_negative_validate_instances(cls: ChemicalClass, dataset: Dataset) -> Tuple[List[ChemicalStructure], List[ChemicalStructure]]:
    """
    Get positive and negative instances for a chemical class from validation set.

    Args:
        cls:
        dataset:

    Returns:

    """
    s2i = dataset.smiles_to_instance()
    all_validation = set(dataset.validation_examples)
    all_positive = set(cls.all_positive_examples)
    all_smiles = dataset.all_smiles()
    positive_examples = list(all_positive.intersection(all_validation))
    positive_instances = [s2i[smiles] for smiles in positive_examples]
    negative_examples = list((all_smiles - all_positive).intersection(all_validation))
    negative_instances = [s2i[smiles] for smiles in negative_examples]
    return positive_instances, negative_instances

def evaluate_class(rset: ResultSet, dataset: Dataset, threshold: float=None) -> EvaluationResult:
    """
    Evaluate a chemical class.

    Args:
        rset: ResultSet containing the best result for a chemical class
        dataset:

    Returns:

    """
    br = rset.best_result
    cls_lite = br.chemical_class
    cls_name = cls_lite.name
    cls = dataset.get_chemical_class_by_name(cls_name)
    logger.info(f"Evaluating {cls.name}")
    positive_instances, negative_instances = get_positive_and_negative_validate_instances(cls, dataset)
    logger.info(f"Validate POS={len(positive_instances)} NEG={len(negative_instances)}")
    code_str = br.code
    test_result = evaluate_program(code_str, cls_lite, positive_instances, negative_instances,
                                   threshold=threshold,
                                   model=br,
                                   )
    logger.info(f"Test Result F1: {test_result.f1} Train F1: {br.f1}")
    return EvaluationResult(train_results=rset, test_result=test_result)


def learn_ontology_iter(dataset: Dataset, config: Config, working_dir: Optional[Path] = None, include_only: Optional[List[str]] = None, exclude: Optional[List[str]] = None, mapped_only = False) -> Iterator[ResultSet]:
    """
    Learn a set of chemical classes.

    Args:
        dataset:
        config:

    Returns:

    """
    expt = config.experiment_name
    cache_dir = working_dir / expt / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    for cls in dataset.classes:
        logger.info(f"Learning {cls.name}")
        safe_cls_name = safe_name(cls.name)
        if include_only:
            if cls.name not in include_only and cls.id not in include_only and safe_cls_name not in include_only:
                logger.info(f"Skipping {cls.name}, not in include_only")
                continue
        if exclude:
            if cls.name in exclude or cls.id in exclude or safe_cls_name in exclude:
                logger.info(f"Skipping {cls.name}, in exclude")
                continue
        if mapped_only:
            if not cls.xrefs:
                logger.info(f"Skipping {cls.name}, not mapped")
                continue
        filename = cache_dir / f"{safe_cls_name}.json"
        if filename.exists():
            logger.info(f"Skipping {cls.name}")
            with open(filename, "r") as f:
                rset = ResultSet.model_validate_json(f.read())
                # yield rset
        else:
            pos, neg = get_positive_and_negative_train_instances(cls, dataset)
            if len(pos) < config.min_positive_examples_for_training:
                logger.info(f"Skipping {cls.name}, {len(pos)} is not enough positive examples")
                continue
            if len(neg) < config.min_negative_examples_for_training:
                logger.info(f"Skipping {cls.name}, {len(neg)} is not enough negative examples")
                continue
            rset = learn_program(cls, pos, neg, config)
            with open(filename, "w") as f:
                f.write(rset.model_dump_json(indent=2))
        yield rset


def apply_threshold(rset: ResultSet, threshold: float) -> ResultSet:
    """
    Apply a threshold to a ResultSet.

    Args:
        rset: ResultSet to apply the threshold to
        threshold: threshold to apply

    Returns:
        ResultSet with the threshold applied
    """
    new_results = []
    for result in rset.results:

        if result.f1 is not None and result.f1 >= threshold:
            new_results.append(result)
    return ResultSet(results=new_results, chemical_class=rset.chemical_class, config=rset.config)