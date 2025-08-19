import math
import re
from copy import copy

from pydantic import BaseModel, Field
from typing import List, Optional, Tuple, Union, Set, Dict

SMILES_STRING = str


# make this hashable?
class ChemicalStructure(BaseModel):
    """Represents a chemical entity with a known specific structure/formula."""
    name: str = Field(..., description="rdfs:label of the structure in CHEBI")
    smiles: SMILES_STRING = Field(..., description="SMILES string derived from CHEBI")

    def __hash__(self):
        return hash(self.smiles)

    def __eq__(self, other):
        if not isinstance(other, ChemicalStructure):
            return NotImplemented
        return self.smiles == other.smiles


class ChemicalClass(BaseModel):
    """Represents a class/grouping of chemical entities."""
    id: str = Field(..., description="id/curie of the CHEBI class")
    name: str = Field(..., description="rdfs:label of the class in CHEBI")
    definition: Optional[str] = Field(None, description="definition of the structure from CHEBI")
    parents: Optional[List[str]] = Field(default=None, description="parent classes")
    xrefs: Optional[List[str]] = Field(default=None, description="mappings")
    all_positive_examples: List[SMILES_STRING] = []

    def lite_copy(self) -> "ChemicalClass":
        """
        Create a copy of the chemical class without the instance fields
        Returns:
        """
        cc = copy(self)
        cc.all_positive_examples = []
        #cc.train_positive = []
        #cc.train_negative = []
        #cc.validate_positive = []
        #cc.validate_negative = []
        return cc

class Dataset(BaseModel):
    """
    Represents a dataset of chemical classes.
    """
    ontology_version: Optional[str] = None
    min_members: Optional[int] = None
    max_members: Optional[int] = None
    classes: List[ChemicalClass]
    structures: List[ChemicalStructure] = None
    validation_examples: Optional[List[SMILES_STRING]] = None

    @property
    def name(self):
        return f"bench-{self.ontology_version}-{self.min_members}-{self.max_members}"

    def all_smiles(self) -> Set[SMILES_STRING]:
        return {s.smiles for s in self.structures}

    def smiles_to_instance(self) -> Dict[SMILES_STRING, ChemicalStructure]:
        return {s.smiles: s for s in self.structures}

    def get_chemical_class_by_id(self, class_id: str) -> ChemicalClass:
        for cc in self.classes:
            if cc.id == class_id:
                return cc
        raise ValueError(f"Class {class_id} not found in dataset")

    def get_chemical_class_by_name(self, class_name: str) -> ChemicalClass:
        for cc in self.classes:
            if cc.name == class_name:
                return cc
        raise ValueError(f"Class {class_name} not found in dataset")

from typing import Optional


class Config(BaseModel):
    """Experimental setup"""
    llm_model_name: str = "gpt-4o"
    experiment_local_name: Optional[str] = None
    f1_threshold: float = 0.8
    max_attempts: int = 4
    min_positive_examples_for_training: int = 23
    min_negative_examples_for_training: int = 23
    min_positive_examples_for_validation: int = 5
    min_negative_examples_for_validation: int = 5
    max_positive_instances: Optional[int] = None
    max_positive_to_test: Optional[int] = None
    max_negative_to_test: Optional[int] = None
    max_positive_in_prompt: int = 50
    max_negative_in_prompt: int = 20
    max_examples_in_feedback: Optional[int] = 25
    test_proportion: float = 0.2
    use_definitions: bool = True
    use_the_force: bool = False

    @property
    def experiment_name(self):
        ln = self.experiment_local_name or "undef"
        model_name = self.llm_model_name
        if "/" in model_name:
            model_name = model_name.split("/")[-1]
        return f"{model_name}-{ln}"


OUTCOME = Union[Tuple[str, Optional[str]], Tuple[str, str, Optional[str]]]

class Outcome(BaseModel):
    smiles: SMILES_STRING
    name: Optional[str] = None
    reason: Optional[str] = None

    class Config:
        frozen = True  # Makes the model immutable

    def __str__(self):
        return f"{self.smiles} ({self.name}): {self.reason}"

    def __repr__(self):
        return f"{self.smiles} ({self.name}): {self.reason}"

    def __hash__(self):
        return hash((self.smiles, self.name, self.reason))


class CodeStatistics(BaseModel):
    """Code statistics"""
    lines_of_code: int
    log_lines_of_code: float
    indent_by_line: List[int]
    max_indent: int
    imports: List[str]
    imports_count: int
    methods_called: List[str]
    methods_called_count: int
    smarts_strings: List[str]
    smarts_strings_count: int
    defs: List[str]
    defs_count: int
    returns: List[str]
    returns_count: int
    complexity: float

    @classmethod
    def from_code(cls, code: str):
        """Extract statistics from code"""
        imports = []
        lines = []
        ignored_lines = []
        indent_by_line = []
        returns = []
        defs = []
        last_line_indent = 0
        in_def = False
        for line in code.split("\n"):
            if line.startswith("__metadata__"):
                break
            num_spaces = len(line) - len(line.lstrip())
            line_indent = num_spaces // 4
            if num_spaces % 4:
                # likely not a "true" new line
                line_indent = last_line_indent
            elif line_indent > last_line_indent + 1:
                # likely a continuation of the previous line
                line_indent = last_line_indent
            last_line_indent = line_indent
            if in_def:
                if line.strip() and line_indent == 0:
                    pass
                    #in_def = False
                else:
                    if line.strip():
                        lines.append(line)
                    indent_by_line.append(line_indent)
            elif line.startswith("from") or line.startswith("import"):
                imports.append(line)
            elif line.startswith("def"):
                in_def = True
                defs.append(line.replace("def", "").strip())
            else:
                ignored_lines.append(line)
        lines_of_code = len(lines)
        max_indent = max(indent_by_line) if indent_by_line else 0
        imports_count = len(imports)
        method_re = re.compile(r"\.(\w+)\(")
        # TODO: this misses when a variable is passed
        smarts_re = re.compile(r"Chem.MolFromSmarts\((.+)\)")
        methods_called = []
        smarts_strings = []
        def de_quote(s):
            if s.startswith("'") and s.endswith("'"):
                return s[1:-1]
            if s.startswith('"') and s.endswith('"'):
                return s[1:-1]
            return s
        for line in lines:
            methods_called.extend(method_re.findall(line))
            smarts_strings.extend([de_quote(x) for x in smarts_re.findall(line)])
            line_lstrip = line.lstrip()
            if line_lstrip.startswith("return"):
                returns.append(line.replace("return", "").strip())
            if line_lstrip.startswith("def"):
                defs.append(line.replace("def", "").strip())
        methods_called = list(set(methods_called))
        methods_called_count = len(methods_called)
        smarts_strings = list(set(smarts_strings))
        smarts_strings_count = len(smarts_strings)
        defs_count = len(defs)
        returns_count = len(returns)
        log_loc = math.log(lines_of_code) if lines_of_code else 0
        complexity = (methods_called_count + defs_count + returns_count + log_loc + max_indent) / 5
        return cls(
            lines_of_code=lines_of_code,
            log_lines_of_code=log_loc,
            indent_by_line=indent_by_line,
            max_indent=max_indent,
            imports=imports,
            imports_count=imports_count,
            methods_called=methods_called,
            methods_called_count=methods_called_count,
            smarts_strings=smarts_strings,
            smarts_strings_count=smarts_strings_count,
            defs=defs,
            defs_count=defs_count,
            returns=returns,
            returns_count=returns_count,
            complexity=complexity,
        )


class Result(BaseModel):
    """Result of running workflow on a chemical class"""
    chemical_class: ChemicalClass
    config: Optional[Config] = None
    code: str
    code_statistics: Optional[CodeStatistics] = None
    message: Optional[str] = None
    true_positives: Optional[List[Outcome]] = None
    false_positives: Optional[List[Outcome]] = None
    true_negatives: Optional[List[Outcome]] = None
    false_negatives: Optional[List[Outcome]] = None
    sample_true_negatives: Optional[List[Outcome]] = None
    sample_false_negatives: Optional[List[Outcome]] = None
    attempt: int = 0
    reasoning: Optional[str] = None
    success: bool = True  ## True if no runtime errors or compilation errors
    best: bool = False
    error: Optional[str] = None
    stdout: Optional[str] = None

    num_true_positives: Optional[int] = None
    num_false_positives: Optional[int] = None
    num_true_negatives: Optional[int] = None
    num_false_negatives: Optional[int] = None
    num_negatives: Optional[int] = None

    precision: Optional[float] = None
    recall: Optional[float] = None
    f1: Optional[float] = None
    accuracy: Optional[float] = None
    negative_predictive_value: Optional[float] = None

    def calculate(self):
        """Calculate derived statistics"""
        self.num_true_positives = len(self.true_positives or [])
        self.num_false_positives = len(self.false_positives or [])
        if self.num_true_negatives is None:
            self.num_true_negatives = len(self.true_negatives or [])
        if self.num_false_negatives is None:
            self.num_false_negatives = len(self.false_negatives or [])
        if self.num_true_positives + self.num_false_positives:
            self.precision = self.num_true_positives / (self.num_true_positives + self.num_false_positives)
        else:
            self.precision = 0.0
        if self.num_true_positives + self.num_false_negatives:
            self.recall = self.num_true_positives / (self.num_true_positives + self.num_false_negatives)
        else:
            self.recall = 0
        if self.precision and self.recall:
            self.f1 = 2 * (self.precision * self.recall) / (self.precision + self.recall)
        else:
            self.f1 = 0
        if self.num_true_positives + self.num_true_negatives + self.num_false_positives + self.num_false_negatives:
            self.accuracy = (self.num_true_positives + self.num_true_negatives) / (
                    self.num_true_positives + self.num_true_negatives + self.num_false_positives + self.num_false_negatives)
        if self.num_true_negatives + self.num_false_negatives:
            self.negative_predictive_value = self.num_true_negatives / (self.num_true_negatives + self.num_false_negatives)
        else:
            self.negative_predictive_value = 0
        if self.code and not self.code_statistics:
            self.code_statistics = CodeStatistics.from_code(self.code)

class ResultSet(BaseModel):
    """A set of results for a chemical class"""
    best_result: Optional[Result] = None
    results: List[Result]
    sorted_attempts: List[int] = []
    experiment_name: Optional[str] = None

    @classmethod
    def from_results(cls, results: List[Result]) -> "ResultSet":
        """Populate the result set from a list of results"""
        obj = cls(results=results, best_result=max(results, key=lambda r: r.f1))
        if obj.best_result:
            obj.best_result.best = True
        return obj

    @property
    def markdown(self):
        """Generate markdown for the result set"""
        br = self.best_result
        best_attempt = br.attempt if br else -1
        chem = br.chemical_class
        md = f"# Results for {chem.id} {chem.name}\n\n"
        for r in self.results:
            md += f"## Attempt {r.attempt}\n\n"
            if r.attempt == best_attempt:
                md += f"**Best result**\n\n"
            if r.message:
                md += f"### Feedback from previous attempt\n\n"
                md += f"{r.message}\n\n"
            if r.reasoning:
                md += f"### Reasoning\n\n"
                md += f"{r.reasoning}\n\n"

            md += "### Code\n\n"
            md += f"```python\n{r.code}\n```\n\n"
            if r.error:
                md += f"### Error\n\n"
                md += f"```python\n{r.error}\n```\n\n"
            md += f"Precision: {r.precision:.2f}\n\n"
            md += f"Recall: {r.recall:.2f}\n\n"
            md += f"F1: {r.f1:.2f}\n\n"
        return md

class EvaluationResult(BaseModel):
    """Result of evaluating a model"""
    train_results: ResultSet
    test_result: Result

    def calculate_reward(self):
        """Calculate derived statistics"""
        tr = self.test_result
        tr.calculate()
        br = self.train_results.best_result
        br.calculate()
        tp_reward = tr.num_true_positives * br.precision
        fp_penalty = tr.num_false_positives * br.precision
        fn_penalty = tr.num_false_negatives * br.negative_predictive_value
        reward = (tp_reward - fp_penalty) - fn_penalty
        return reward

    @property
    def markdown(self):
        """Generate markdown for the evaluation result"""
        return self.train_results.markdown


class EvaluationExperiment(BaseModel):
    """Represents an evaluation experiment"""
    config: Config
    evaluation_results: List[EvaluationResult]

class ClassificationResult(BaseModel):
    input_smiles: SMILES_STRING
    class_id: str
    class_name: Optional[str] = None
    is_match: bool
    reason: Optional[str] = None
    confidence: Optional[float] = None
