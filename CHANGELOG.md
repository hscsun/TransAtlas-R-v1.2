# Changelog

All notable changes to TransAtlas are recorded here. The project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.2.0] — 2026-11-12

### Added

- Preserve every activity profile in differential tables, with sample counts,
  test status and reason. Undefined Welch tests retain descriptive statistics
  and NA P values; BH correction uses only successfully tested profiles.
- Record each run warning with stage, call, class and timestamp in warnings.csv,
  with per-stage records, warning_summary.csv and status_warnings.json.
- Record actual MSigDB database/version, collection, target species, ortholog
  mapping and gene-set membership checksum for GSEA and GSVA.
- Recognize mouse-native Hallmark as MH. H keeps the human Hallmark ortholog
  policy for compatibility; MH explicitly requests mouse-native Hallmark.
- Remove repeated merged-guide overrides and account for intentionally unlabelled
  GO hub nodes. Empty chromosome mappings record NA instead of an empty maximum.


Revision on 2026-11-12: every calculation now uses unbounded log2FC.
`up_cutoff` (2.5) controls gene displays only; `low_cutoff` (0.15, NULL/0
disables) controls GSEA independently of the screening ladder. ORA backgrounds
use the eligible tested universe. Negative duplicate-ID effects are selected
by largest absolute magnitude. TF module target summaries count each target
once per TF. GSEA, ORA and regulatory results are expected to change.
The reported TF odds ratio uses the full 2x2 contingency table; chromosome
thresholds are included exactly even with a custom screening ladder. Missing
network effects retain NA, and reports exclude nonsignificant GO fallback rows.
The standalone GO network uses unshrunk estimates by default and supports the
configured display cutoff; its stress layout uses graphlayouts. The localization
report uses the annotated DEG denominator, and key-target table references and
shrinkage descriptions match the current report.

### Fixed

- Read CSV/Excel/TSV and data.frame metadata through one validated path; automatic
  discovery skips mismatched sample subsets.
- Parse CLI threshold vectors and validate effect-size/significance parameters;
  include the primary effect threshold in every screening ladder.
- Invalidate DESeq2 cache for metadata assignments, alpha, filtering, shrinkage,
  model/contrast settings and package versions. Re-run requested QC on cache hits.
  Recreate upstream CSV exports from the cached fit, including when the output
  folder initially contains only the cache.
- Fall back from vst to varianceStabilizingTransformation for small gene sets.
- Handle missing regulon GSEA columns, missing package data, empty length mappings,
  and GFF3 exon-to-transcript-to-gene relationships.
- Preserve database exon-union FPKM fallback; record source and coverage. An active
  FPKM floor stops clearly when no usable FPKM matrix can be built. NULL disables it.
- Separate significant GO tables from all-term fallback tables when reloading.
- Use the project random seed and project-relative cache paths; load sections in
  isolated environments and validate returned contexts.
- Record optional errors, QC/section/report failures and an overall run status;
  protect graphics devices and diagnose invalid/failed STRING batches.
- Fix missing-image HTML rendering, section status slugs and ranking descriptions.
- Use consistent strict adjusted-P and inclusive effect-size boundaries, shared
  alpha defaults, and explicit exploratory PPI annotation and plot label parameters.

### Changed

- The volcano and expression-versus-effect panels no longer count the genes
  below the expression floor in their legend; the floor class stays labelled.
  `S1_expression_effect` is drawn at `c(6.5, 6.0)` and `s3_msigdb_barplot` at
  `c(12, 8)`.
- The GO term–gene network of `scripts/go_term_gene_network.R` uses the
  Kamada–Kawai layout, labels more genes, carries the ontology term names as
  radial labels with leader lines, and keeps only the log2FC colour bar, placed
  inside the lower-right corner.
- Type sizes raised by one point in several panels (GO term–gene network, TF
  dot plot, TF–target network, TF modules, section 9 heat maps, PPI hub panel),
  and the string `HALLMARK_` is stripped from the MSigDB labels.

### Documentation

- README rewritten in English throughout, with author contributions stated.
- `docs/脚本来源_特色与命名分析.md` retired.

## [1.1.0] — 2026-10-26

One effect-size cut-off and one expression floor, shared by every section.

### Changed

- `log2fc` (default **0.75**) is now the only |log2FC| cut-off. It calls a gene
  differentially expressed and builds the ranked lists and gene sets of the
  enrichment and regulatory sections; `gene_set_lfc_threshold` is retired.
- `min_fpkm` (default **5**) is the expression floor of the whole workflow: the
  higher of a gene's two group-mean FPKM values must reach it, otherwise the
  gene is dropped from the contrast table before screening. It replaces the
  transcription-factor-only `tf_min_fpkm` (10), which is retired, and it is
  applied once in the run context so every section sees the same universe. The
  table behind it is written to
  `00_differential_expression/expression_floor.csv`.
- Both retired arguments are still accepted and produce a deprecation warning;
  the command-line flags `--gene-set-lfc-threshold` and `--tf-min-fpkm` behave
  the same way. Use `--log2fc` and `--min-fpkm`.
- A typical run now needs four settings: `species`, `log2fc`, `min_fpkm` and
  `log2fc_cap`.
- The legacy `log10(baseMean) <= 2` rule of `AnaDGEsS1.R` is off by default
  (`expression_log10_threshold = NULL`). It was a second, count-based
  expression filter that removed short genes reaching the FPKM floor and kept
  long genes below it; `min_fpkm` is now the only expression threshold. Pass
  `expression_log10_threshold = 2` to restore the old behaviour.
- The case study was re-run with the new defaults: 8,369 genes pass the
  expression floor, 1,213 are differentially expressed at |log2FC| ≥ 0.75.

### Added

- `00_differential_expression/expression_floor.csv` records the group-mean FPKM
  maximum of every gene and whether it was kept, and `status.json` files record
  the floor that was applied.

## [1.0.0] — 2026-10-23

First public release, accompanying the software paper *TransAtlas: a
reproducible eleven-section R workflow that links differential expression to
regulatory, interaction and localisation evidence in bulk RNA-seq*.

### Analysis

- A single entry point, `run_all()`, runs eleven paper-ordered sections from a
  count workbook and a two-column sample table: screening with a
  |log2FC| ladder (0.15–2.5), GO over-representation with semantic
  simplification and ten functional classes, threshold-free enrichment (GO and
  MSigDB Hallmark), KEGG with `pathview` maps, transcription-factor regulon
  inference from DoRothEA, regulon activity (decoupleR ULM/AUCell/WMean),
  TF–target networks, TF modules by Jaccard overlap, gene-set heat maps,
  STRING protein–protein interaction modules and subcellular localisation.
- Species presets for mouse, human and rat set the OrgDb, KEGG code,
  identifier types, MSigDB collection and DoRothEA collection together.
- FPKM is computed from the exon-union gene lengths carried by the count
  workbook (`Exonic.gene.sizes` by default), with a GTF/GFF or the species
  annotation package as fallback; a pre-computed expression matrix is used
  when one is supplied.
- Every section writes a `status.json`, machine-readable tables and figures at
  fixed physical sizes, so a single section can be re-run on its own.
- A provenance record captures package versions, the platform and the exact
  resolved call; one run also produces a figure atlas and a fixed-template
  Word report.

### Case study

- Serum restimulation of NIH 3T3 fibroblasts (10% versus 2% FBS, three
  libraries per arm): 55,401 quantified genes, 14,655 tested, 807 DEGs at
  adjusted P < 0.05 and |log2FC| ≥ 1, 1,509 at |log2FC| ≥ 0.75.
- The case-study inputs (`input/NIH3T3.xlsx`, `input/NIH3T3_metadata.csv`) and
  the scripts that assemble the manuscript (`work/`) are part of the release.

### Documentation

- `README.md` covers installation, the input contract, the parameter set, the
  output layout and citation.
- `docs/USER_GUIDE.md` walks through one complete run with the NIH 3T3 example:
  setup, the run command, how to verify the result, which file answers which
  question, how to change the analysis and what to do when a section fails.
- `CITATION.cff` and `.zenodo.json` provide the citable record.
