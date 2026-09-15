# RTM3

This simulation study examines how seizure-detector sensitivity and false-alarm rate affect trial eligibility, observed regression to the mean, and apparent placebo response.

The main analysis uses a two-month baseline, a three-month test period, and 30-day months. Eligibility requires at least four seizures/month on average during baseline and at least three in each baseline month. It omits the former 25-day seizure-free-interval rule. No treatment effect is imposed.

## Reproduce the paper and appendix

Use Python 3.12 in a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/reproduce_all.py
```

The last command executes both notebooks. The main-paper notebook calls the canonical analysis, generating 100,000 potential participants, then creates both manuscript figures through shared functions. The Appendix notebook reads that run's aggregate results and embeds all sensitivity analyses, equations, tables, and figures. HTML copies are written to `results/notebooks/` for review without Jupyter.

Use `requirements-reproduction.txt` in place of `requirements.txt` to pin the direct dependencies to the released run. HTML equations use MathJax and may need network access; notebook equations render in Jupyter.

For the analysis alone:

```bash
python reproduce_rtm3.py
python -m unittest discover -s tests -v
```

For a short test run that does not overwrite the paper results:

```bash
python reproduce_rtm3.py --n-patients 1000 --output-dir results_smoke --progress-every 0
```

## Sources and outputs

- `reproduce_rtm3.py`: canonical simulation, eligibility, endpoints, diagnostics, validation, and Figure 1.
- `realSim.py`: unchanged copy of the CHOCOLATE simulator source. Provenance and SHA-256 are recorded in the run metadata.
- `rtm3_reporting.py`: notebook tables, comparison figures, interpretation, and conceptual Figure 2.
- `for_rtm3_paper.ipynb`: executed full reproduction with both manuscript images embedded and the fixed-versus-reselected cohort comparison quoted in the main text.
- `for_appendix_RTM3.ipynb`: executed methodological appendix with every condition, including zero-baseline denominators.
- `appendix_methods.md`: equations and methods used to build the appendix.
- `scripts/build_paper_notebook.py` and `scripts/build_appendix_notebook.py`: notebook builders; pass `--execute` to either.
- `tests/test_reproduce_rtm3.py`: regression tests for shared participant traits, strict RTM ties, reference scales, eligibility, and reproducibility.

The `results/` directory contains primary, interval, correction, selection, fixed-cohort, and additive-threshold CSVs. It also contains PNG/TIFF figures, `rtm3_validation.json`, and `rtm3_run_metadata.json`. Metadata record seeds, software versions, source hashes, output hashes, and simulator provenance. Figure 2 is a conceptual diagram, not an estimated causal model.

The main seed is `20260813`. Detection and false alarms use `20260814` and `20260815`. The five-month trial diary and separate 36-month reference realization share sampled rate, cluster status, cycle frequencies, and amplitudes. Detector conditions share latent participants; detections and false alarms are nested within each sweep. Exact scaled integer comparisons prevent equal-distance RTM ties from being misclassified. Package versions used for the released run are recorded in metadata; different package versions can change pseudorandom results or image rendering.

## Interpretation and sensitivity analyses

The observed-RTM indicator requires baseline above the detector-scale long-term reference and a test period strictly closer to that reference. It does not identify an isolated causal "RTM Type 3" component.

The revised 100,000-participant run gives perfect-detector placebo median percentage change of 14.3%, increasing to 48.1% at 10% sensitivity and 31.2% at one false alarm/day when eligibility is re-evaluated. The rising observed-RTM pattern reverses when the perfect-detector cohort is held fixed. Thus reselection matters to the primary result. Fixed-cohort percentage change is undefined at zero measured baseline; the appendix reports those participants and the valid denominator rather than assigning a pseudocount.

Using raw-count eligibility thresholds with expected alarms added selects exactly the same participants as subtracting expected alarms from counts. The separate monthly thresholds are `3 + 30*FAR`, and the baseline mean threshold is `4 + 30*FAR`. The two-month total equivalent is `8 + 60*FAR`. Outcomes remain corrected in both analyses. This is an equivalent rule, not a more precise false-alarm estimator.

The uncorrected-outcome comparison is different: it uses raw counts and the reference `sensitivity * effective_long_term_frequency + 30*FAR`. Exact oracle removal uses known simulated alarm labels and is not an available correction for unlabeled clinical events.

The original 25-day rule and alternative interval implementations are appendix sensitivities. Daily sampling cannot represent 12-hour or 24-hour within-day cycles. The sweeps vary sensitivity and FAR separately and do not predict outcomes for a specific device at paired operating characteristics.

## Historical files

`old_for_rtm3_paper.ipynb` preserves the previous main notebook. Older notebooks, `rtm_tools*` helpers, and older result files are historical, not canonical reproduction sources. Their computations and figures should not be mixed with this revision.

In the repository, the Word-verified submission files are in `submission/epileptic_disorders/`, with a standalone `Reproducibility/` package. DOCX files preserve live Zotero citations. Citation refreshes, final saves, and page checks were completed in Microsoft Word and Zotero; future document edits require those checks again. The document-revision helpers are archival transformations, separate from the canonical simulation workflow.

## License

See `LICENSE`.
