# Revision status

The core simulation, two executed notebooks, tables, and figures are complete. All 12 regression tests pass. A standalone reproduction in a clean temporary directory produced byte-identical CSVs and manuscript PNG/TIFF files. Notebook images are embedded, and the main notebook's PNG bytes match the generated files. The DOCX images match the regenerated assets.

Validated code, notebooks, and aggregate outputs were committed and pushed to origin/main as `49cdcc62f1d5451a43d6f1de3074e9a10db3f28d`. The local and remote commit hashes match. Unverified DOCX drafts and unrelated pre-existing changes were not included in that commit.

The scientific conclusion needs the new qualification now included in the draft. Reselected observed RTM increases from 37.2% with perfect detection to 64.4% at 10% sensitivity and 53.6% at one false alarm/day. Holding enrollment fixed gives 29.3% and 34.7%, respectively. Fixed-cohort MPC uses positive-baseline participants and excludes 6,182 and 538 undefined values at those two conditions. Additive raw-count thresholds are exactly equivalent to corrected-count eligibility. Retaining the original 25-day rule changes MPC by at most 2.29 percentage points.

The unchanged upstream simulator was verified against the local CHOCOLATE file, with its commit and SHA-256 recorded in metadata. The previous main notebook was preserved as `old_for_rtm3_paper.ipynb`. The new notebook calls shared functions and reruns the full 100,000-participant analysis. Direct dependency versions are pinned in `requirements-reproduction.txt`.

The remaining Word work is complete. All four DOCX files were saved and inspected in Microsoft Word. The final 13-page manuscript PDF was exported by Word and every page inspected. The abstract and captions now bold only their labels, the key points are boxed, the educational questions are kept together, and the title page fits on one page. The standalone disclosure matches the manuscript with abbreviations expanded.

Zotero refreshed the live citation fields and bibliography. The document uses stored journal abbreviations rather than automatic MEDLINE abbreviation of extended journal titles. The missing journal abbreviation in the older Zaccara item AX5H9B4Q was set to `Epilepsy Behav` in Zotero. The final bibliography correctly abbreviates references 1, 9, and 18. Ordered citation identities, the two embedded figures, and all nonbibliographic manuscript prose survived the refresh unchanged. No tracked insertions or deletions remain.

The main-paper notebook now displays the fixed-versus-reselected conditions quoted in the manuscript. Both notebooks were rerun. A final independent audit confirmed the fixed-cohort numbers and causal qualifications. It found one overgeneralized Appendix sentence about alternative interval implementations. This was corrected: observed-day intervals change RTM and MPC less than the original rounded-count interval rule, while oracle-day intervals change them more, comparing maximum absolute differences across the sweep.

The original iCloud submission documents were checked against source backups before replacement. The prior set is retained in `Archive_20260915`; the current four DOCX files, manuscript PDF, two TIFFs, and complete Reproducibility folder are the revised set. No journal submission was performed. Prior scientific-meeting presentation remains an author-supplied fact; no statement was invented. Unrelated pre-existing deletions and edits were left untouched.

The notebook kernel printed a sandbox permission warning during shutdown when checking child processes. Execution returned success, every code cell has an execution count, no cell contains an error output, and all saved results passed validation. The warning did not affect calculations or figures.

The `unslop` writing check removed ambiguous causal claims, undefined initials, and redundant phrasing. A separate manuscript-to-results audit found no numerical discrepancy; its final two copyedits were applied.
