# PAT-style review: Epilepsia Brief Communication draft

Manuscript reviewed: `epilepsia_brief_communication_draft.md`  
Review date: May 28, 2026  
Target journal/type: Epilepsia Brief Communication

## Deterministic metrics

Post-revision `pat_metrics.py` output:

- Total manuscript words including title page, summary, declarations, references, and legends: 2,446.
- Summary: 180 words.
- Main text: 1,299 words.
- References: 14.
- Figures/tables: two figures, zero tables.
- Average sentence length: 16.6 words.
- Sentences >30 words: 10.9%.
- Flesch-Kincaid grade: 14.4.
- Passive voice estimate: 23.1%.

The long-sentence warning is mostly caused by Markdown superscript citation tags preventing sentence splitting in the metrics script. The scientific prose is within the expected academic range for a short methods-heavy manuscript.

## Journal compliance audit

Pass:

- Brief Communication length: main text is below the 1,800-word guideline.
- Summary: unstructured and below the 200-word Brief Communication limit.
- Short summary: provided in the manuscript and separately as `epilepsia_brief_communication_short_summary.md`.
- Keywords: six terms provided.
- References: 14, below the limit of 15, and ordered numerically.
- Figures/tables: two total display items, within the Brief Communication limit.
- Figure format: submission TIFFs created as `Goldenholz-fig1.tif` and `Goldenholz-fig2.tif`, RGB, LZW compressed, 600 dpi metadata, below 40 megapixels.
- Figure legends: both legends explain symbols/abbreviations and identify source PNG files.
- Word draft: `epilepsia_brief_communication_draft.docx` generated with 1-inch margins, Times New Roman 12-point body text, double spacing, and first-author/page-number header.
- Funding, disclosure, data/code availability, author contributions, and ethical publication statement are present.
- Human/animal subjects language is appropriate for a synthetic simulation study.

Items requiring author confirmation before submission:

- Corresponding author address, telephone, fax, and email are placeholders.
- Author affiliations and author contribution text should be confirmed.
- RB and BW conflict-of-interest statements are drafted as no conflicts but require confirmation.
- Text page count should be finalized after Word formatting.
- References should be imported into a reference manager or checked against PubMed/Wiley before submission, although the main high-priority references were spot-checked and the CHOCOLATES author line now matches PubMed.

## Scientific/reviewer-risk audit

Main strengths:

- The manuscript has a clear single-message brief communication: device sensitivity and false alarm rate can alter RTM and placebo MPC even before any active treatment effect.
- The methods now explicitly restrict eligibility testing to baseline only, matching the revised analysis scope.
- The figure legend and Methods state that the error bars are 95% confidence intervals for the estimated median MPC, not confidence regions for individual percentage change.
- Limitations address simulation abstraction, uniform detector performance, idealized FAR correction, and CI interpretation.

Likely reviewer questions:

- Whether daily independent sensitivity and Poisson FAR are too simple for real detector behavior.
- Whether focal seizure detection performance should be distinguished from tonic-clonic wearable-device performance.
- Whether the idealized expected-FAR subtraction is clinically implementable.
- Whether device nonwear, missingness, and seizure-type-specific sensitivity should be simulated in future work.

No blocking scientific issue was identified for a draft submission. The remaining blockers are administrative placeholders and final reference-manager verification.

## Revisions made after review

- Checked the CHOCOLATES citation against PubMed and restored the two-author listing.
- Replaced preliminary figure TIFF names with Epilepsia-style first-author figure filenames.
- Recreated figure TIFFs at 600 dpi with LZW compression.
- Added a separate short-summary supporting document.
- Updated the manuscript word-count line.
- Generated and rendered a Word draft; visual QA found no obvious clipping, overlap, or missing text.
