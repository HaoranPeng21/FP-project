# Local verification

Verified on 2026-11-03 using the existing private project inputs and the R
environment in `DEPENDENCIES.md`.

**Ready within the reviewed data-to-figure scope:** 19 of 19 analysis scripts
completed successfully and generated 27 PDF panels/panel groups. All 27 R files
parsed successfully. Missing-data behaviour is explicit; no example or synthetic
records were substituted for the private inputs.

## Numerical reconciliation

| Check | Result |
|---|---|
| Full-catalogue Sylph benchmark | TP 4,453; FN 47; FP 1; TN 70,282 |
| Shared24 Sylph benchmark | TP 1,293; FN 7; FP 0; TN 20,324 |
| Shared24 MetaPhlAn4 benchmark | TP 1,299; FN 1; FP 849; TN 19,475 |
| SGB15332 | TP 63; FN 0; FP 581; TN 257 |
| Pathway tests | All 5,256 Fisher tests/BH values reconciled; all original star-eligibility decisions retained |
| Pathway display | 24 selected species–pathway entries; 14 nonredundant profiles plus 2 backgrounds |
| Host composition | 15,811 genomes; percentages sum to 100 within each host |
| Infant summaries | Supplied 10-species × 7-timepoint curves; detected-sample abundance endpoint and frozen 20-test two-part GEE table |
| ABO pooled models | All 14 effect/SE/CI/P/q rows exactly equal the frozen reference table |
| ABO three-comparator models | All 42 effect/SE/CI/P/q rows exactly equal the frozen reference table |
| ABO pairwise tests | All 84 Wilcoxon P/q rows exactly equal the frozen reference table |
| Gene sharing | 25 representatives; 300 pairs; symmetric Jaccard matrix with unit diagonal; rho = −0.72655336534493 |
| Within-species gene display | Frozen PC5 inputs preserve the 8 BMI, 8 tricyclic-antidepressant and 10 PPI families |
| Maast geography | Analyzed n = 1,251 / 631 / 503 / 203; PERMANOVA, PERMDISP and Fitch statistics reconcile to the frozen reference within 10⁻¹⁰; P/q unchanged |

The ABO and Maast comparisons use the existing reference result tables separately
from the new computation. The R scripts also enforce dataset dimensions,
join-key uniqueness and relevant numerical identities at runtime.

## Source discrepancies retained explicitly

The supplied Figure S2D–E titles say **same21**, but their confusion counts and
Figure 1C metrics correspond to **shared24**, with 21,624 cells and 1,300 targets
per tool. The code uses shared24 and corrects its output titles. Same21 inputs are
not substituted for these figure data.

Figure 3A contains 14 species, whereas the supplied 3B–C follows the earlier
ten-species two-part analysis and C shows mean abundance among detected samples.
The later 14-species zero-inclusive trajectory revision is a different analysis
version. The code follows the supplied PDF's panel-specific versions and retains
the original ten-species endpoint-wise GEE correction families for B–C.

The original QAP figure label is P = 10⁻⁵. The native R recomputation with 99,999
permutations gives P = 2×10⁻⁵. The statistic and species-label permutation design
are the same; R and the original NumPy calculation use different random-number
generators. Both P-values are exported, and the original verified value is used
in the figure label. No seed was selected to force agreement.

Figure S1 is an illustrative R reconstruction because its original drawing code
was not found. Its timeline uses public labels from the supplied PDF; its release
crosswalk uses the modal label tuple for each displayed lineage in the private
catalogue. It is not a reconstruction of an unpublished source script.

## Delivery and limits

The release directory passes the code/documentation allowlist and contains no
input data, sample lists, figures, symlinks, private local configuration, machine
home/scratch paths or external executable calls. The private execution receipts,
input checksums, numerical reconciliation and rendered outputs are retained
outside this directory. `.gitignore` excludes data and outputs by default.

Rendering was checked for visible data marks, axes, labels and legends. The code
preserves scientific content and panel groupings; fonts, spacing and manual
manuscript assembly are not guaranteed to match the supplied PDFs pixel for pixel.

These checks validate the stated entry points. They do not rerun upstream read
processing, genome reconstruction, alignment, profiling, pangenome construction,
pyseer fits, longitudinal models or phenotype meta-analysis estimation. Those
scientific products are explicit private inputs rather than hidden script calls.
