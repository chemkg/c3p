RUN = poetry run

test: doctest pytest

pytest:
	poetry run pytest tests

DOCTEST_DIR = c3p
doctest:
	find $(DOCTEST_DIR) -not -path "c3p/programs/*.py" -type f \( -name "*.rst" -o -name "*.md" -o -name "*.py" \) -print0 | xargs -0 $(RUN) python -m doctest --option ELLIPSIS --option NORMALIZE_WHITESPACE


%-doctest: %
	$(RUN) python -m doctest --option ELLIPSIS --option NORMALIZE_WHITESPACE $<

# Experiments

DATASET = results/2025/benchmark/dataset.json
VERBOSE = 1
EX_CORE = claude-3-sonnet-undef o1-undef claude-3-sonnet-force  gpt-4o-undef gpt-4o-hi gemini-2.0-flash-exp-undef o3-mini-undef o3-mini-hi DeepSeek-R1-undef ensemble-9

LEARN = $(RUN) c3p-learn --verbose $(VERBOSE)

2025-llama-405b:
	$(LEARN)  -m lbl/llama-3-405b --dataset $(DATASET) -w results/2025 -x -

2025-llama:
	$(LEARN)  -m lbl/llama --dataset $(DATASET) -w results/2025 -x -

2025-claude-sonnet:
	$(RUN) c3p-learn  --verbose $(VERBOSE)  -m claude-3-sonnet --dataset $(DATASET) -w results/2025 -x -

2025-claude-sonnet-use-the-force:
	$(LEARN) -e force  -m claude-3-sonnet --use-the-force --dataset $(DATASET) -w results/2025 -x -

2025-gpt-4o:
	$(LEARN) -m gpt-4o --dataset $(DATASET) -w results/2025 -x -

2025-o3-mini:
	$(LEARN) -m o3-mini --dataset $(DATASET) -w results/2025 -x -

2025-gemini2-flash:
	$(LEARN) -m gemini-2.0-flash-exp --dataset $(DATASET) -w results/2025 -x -

2025-deepseek-coder:
	$(LEARN) -m deepseek-coder --randomize-order --dataset $(DATASET) -w results/2025 -x -

2025-deepseek-coder-use-the-force:
	$(LEARN) -e force  -m deepseek-coder --use-the-force --dataset $(DATASET) -w results/2025 -x -

2025-cborg-deepthought:
	$(LEARN) -m cborg-deepthought --dataset $(DATASET) -w results/2025 -x -

2025-o1:
	$(LEARN) -m o1 --dataset $(DATASET) -w results/2025 -x -

2025-gpt-4o-hi:
	$(LEARN)  -e hi -m gpt-4o -f 0.9 -a 6 --dataset $(DATASET) -w results/2025 -x -

2025-o3-mini-hi:
	$(LEARN)  -e hi -m o3-mini -f 0.9 -a 6 --dataset $(DATASET) -w results/2025 -x -

2025-deepseek-hi:
	$(LEARN)  -e hi -m deepseek-coder -f 0.9 -a 6 --dataset $(DATASET) -w results/2025 -x -X "carbamate ester" -

2025-deepseek-reasoner:
	$(LEARN)  -m deepseek-reasoner --dataset $(DATASET) -w results/2025 -x  -X "carbamate ester" -X "peptide antibiotic" -X "triglyceride" -

2025-together-r1:
	$(LEARN)  -m deepseek-ai/DeepSeek-R1 --dataset $(DATASET) -w results/2025 -x -

2025-deepseek-reasoner-use-the-force:
	$(LEARN) -e force  -m deepseek-reasoner --use-the-force --dataset $(DATASET) -w results/2025 -x  -X "carbamate ester" -

2025-gpt-4o-nodefs:
	$(RUN) c3p-learn -e nodefs -m gpt-4o --exclude-definitions --dataset $(DATASET) -w results/2025 -x -

2025-r1-distill:
	$(LEARN)  -m groq/deepseek-r1-distill-llama-70b --dataset $(DATASET) -w results/2025 -x -


# Create ensemble model

create-ensemble:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensemble-9 -m 9

create-ensemble-6:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensemble-6 -m 6

create-ensemble-7:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensemble-7 -m 7

create-ensemble-3:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensemble-3 -m 3

create-ensemble-4:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensemble-4 -m 4

eval:
	$(RUN) c3p-validate --verbose 1 -w results/2025/ensembl-5/cache -d $(DATASET)

all-validate: $(patsubst %, validate-%, $(EX_CORE))
validate-%:
	$(RUN) c3p-validate --verbose 1 -w results/2025/$*/cache -d $(DATASET)

summarize-%:
	$(RUN) c3p-summarize --verbose 1 -w results/2025/$*/eval -o results/2025/$*/summary

EX_ORIGINAL = claude-3-sonnet-undef  gpt-4o-undef o1-undef gpt-4o-nodefs gpt-4o-hi ensembl-5 deepseek-coder-undef
EX_COMPLETED = $(EX_ORIGINAL) 
PARTIAL =  gpt-4o-hi  deepseek-coder-undef deepseek-reasoner-undef deepseek-reasoner-force deepseek-chat-hi ensembl-7

compare:
	$(RUN) c3p-compare -x -o results/2025/comparison $(patsubst %,results/2025/%, $(EX_CORE)) results/2025/smartifier

compare-all:
	$(RUN) c3p-compare -x -o results/2025/comparison-partial $(patsubst %,results/2025/%, $(EX_COMPLETED) $(PARTIAL) chebifier)

compare-chebifier:
	$(RUN) c3p-compare -x -o results/2025/comparison-chebifier $(patsubst %,results/2025/%, ensemble-9 chebifier)


compare-all-chebifier:
	$(RUN) c3p-compare -x -o results/2025/comparison-all-chebifier $(patsubst %,results/2025/%, $(EX_CORE) chebifier)


compare-llama:
	$(RUN) c3p-compare -x -o results/2025/comparison-partial results/2025/claude-3-sonnet-undef  results/2025/gpt-4o-undef results/2025/o1-undef results/2025/gpt-4o-nodefs results/2025/gpt-4o-hi results/2025/ensembl-5 results/2025/deepseek-coder-undef results/2025/llama-undef


compare-test:
	$(RUN) c3p-compare -x -o results/2025/comparison-test results/2025/deepseek-reasoner-undef  results/2025/deepseek-reasoner-force

compare-main:
	$(RUN) c3p-compare -x -o results/2025/comparison-main results/2025/ensemble-9 results/2025/chebifier results/2025/smartifier

compare-smartifier:
	$(RUN) c3p-compare -x -o results/2025/comparison-smartifier results/2025/ensemble-9 results/2025/smartifier

compare-reasoners:
	$(RUN) c3p-compare -x -o results/2025/comparison-reasoners results/2025/o1-undef results/2025/o3-mini-undef results/2025/claude-3-sonnet-undef results/2025/deepseek-reasoner-undef 

compare-reasoners2:
	$(RUN) c3p-compare -x -o results/2025/comparison-reasoners2 results/2025/o1-undef results/2025/o3-mini-undef results/2025/DeepSeek-R1-undef 

compare-reasoners3:
	$(RUN) c3p-compare -x -o results/2025/comparison-reasoners3 results/2025/o1-undef results/2025/o3-mini-undef results/2025/DeepSeek-R1-undef results/2025/deepseek-reasoner-undef 

compare-force:
	$(RUN) c3p-compare -x -o results/2025/comparison-force results/2025/claude-3-sonnet-undef results/2025/claude-3-sonnet-force

compare-reasoner-force:
	$(RUN) c3p-compare -x -o results/2025/comparison-reasoner-force results/2025/deepseek-reasoner-undef results/2025/deepseek-reasoner-force

compare-deepseek-hi:
	$(RUN) c3p-compare -x -o results/2025/comparison-reasoner-hi results/2025/deepseek-coder-undef results/2025/deepseek-coder-hi


FIGS = 1 2 3 4 5 6 7 8 9 10 S1 S2 S3 S4 S5
all-figs: $(patsubst %,Figures/Figure-%.jpg,$(FIGS))

sync-figs:
	cp Figures/*jpg Figures-symlink/

Figures/Figure-%.jpg: Figures/Figure-%.png
	convert $< -quality 100 $@

Figures/Figure-1.jpg:
	cp Figures-from-slides/c3po-0.jpg $@

Figures/Figure-2.jpg:
	cp Figures-from-slides/c3po-1.jpg $@

Figures/Figure-3.jpg:
	cp Figures-from-slides/c3po-2.jpg $@

Figures/Figure-4.jpg:
	cp Figures-from-slides/c3po-3.jpg $@

Figures/Figure-5.png:
	cp results/analysis/c3p/individual_c3p/c3p_performance_analysis.png $@

Figures/Figure-6.jpg:
	cp Figures-from-slides/c3po-5.jpg $@

Figures/Figure-7.png:
	magick results/2025/comparison/model_comparison_macro.png results/2025/comparison/model_comparison_micro.png -append $@

Figures/Figure-8.png:
	magick results/2025/comparison-chebifier/model_comparison_macro.png results/2025/comparison-chebifier/model_comparison_micro.png -append $@

Figures/Figure-9.png:
	cp results/2025/comparison-chebifier/method_comparison_ensemble-9_vs_chebifier.png $@

Figures/Figure-10.png:
	cp results/analysis/enrichment_analysis.png $@

Figures/Figure-S1.png:
	cp results/2025/comparison/rank_distribution_top_50.png $@

Figures/Figure-S2.png:
	cp results/2025/ensemble-9/summary/scatterplot-train_f1-f1.png $@

Figures/Figure-S3.jpg:
	cp Figures-from-slides/c3po-13.jpg $@

Figures/Figure-S4.png:
	cp results/2025/comparison/complexity_comparison.png $@

Figures/Figure-S5.jpg:
	cp Figures-from-slides/c3po-15.jpg $@


from-slides:
	convert -density 600 ~/Downloads/C3PO\ Figures.pdf -quality 100 Figures-from-slides/c3po.jpg
