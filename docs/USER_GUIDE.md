# TransAtlas user guide

This guide walks through one complete run with the released example — the
NIH 3T3 serum-restimulation experiment — from installation to reading the
results. Every command below can be copied as it stands; the numbers quoted
are the ones the shipped run produces, so you can tell whether your run
behaved the same way.

If you only want the one-paragraph version, it is:

```r
source("pipeline/run_all.R")
run_all(counts_file = "input/NIH3T3.xlsx",
        samples     = "input/NIH3T3_metadata.csv",
        species     = "mouse",
        contrast    = c("FBS", "Control"))
```

---

## 1. Before you start

| Requirement | Why |
| --- | --- |
| R 4.2 or newer (tested on 4.6.1 / Bioconductor 3.23) | the workflow language and the Bioconductor release it is built on |
| Pandoc 2.11 or newer | renders the HTML manuscript draft and the Word report |
| C/C++ toolchain | a few Bioconductor packages compile from source on first install |
| Network access on the first run | KEGG pathway definitions and STRING interactions are downloaded and then cached |
| About 250 MB of disk per run | figures, tables, the atlas and the report; the shared `cache/` grows to a few GB once KEGG maps and STRING responses are stored |

No GPU is needed. A first run with cold caches takes noticeably longer than a
second one, because sections 4 (KEGG) and 10 (STRING) fetch their resources
once and reuse them from `cache/` afterwards.

## 2. The example dataset

`input/NIH3T3.xlsx` is a single sheet with one row per gene:

| Column | Meaning |
| --- | --- |
| `gene_id` | Ensembl mouse gene identifier (the key the pipeline joins on) |
| `Exonic.gene.sizes` | exon-union length in bp; used to compute FPKM |
| `CON-1`, `CON-2`, `CON-3` | raw counts of the three 2% FBS reference libraries |
| `FBS-1`, `FBS-2`, `FBS-3` | raw counts of the three 10% FBS libraries |
| `GeneSymbol` | mouse gene symbol, used for labels and for GO/KEGG mapping |

`input/NIH3T3_metadata.csv` is the sample table:

```csv
sample_id,condition,replicate
CON-1,Control,1
...
FBS-1,FBS,1
...
```

Three things to notice, because they are what the pipeline reads first:

- the sample columns of the workbook and the `sample_id` values of the table
  must be the **same six strings** — this is checked before anything is fitted;
- `condition` names the two groups; `contrast = c("FBS", "Control")` says which
  way round the comparison runs (FBS over Control, so a positive log2FC means
  higher in 10% FBS);
- the contrast name derived from it, `FBS_vs_Control`, appears in every output
  file name of the run.

## 3. Step 1 — install

```bash
git clone https://github.com/ChangyeSun/TransAtlas-R.git
cd TransAtlas-R
Rscript scripts/install_packages.R
```

The script installs exactly the packages `run_all()` checks at start-up:
Bioconductor packages through `BiocManager`, CRAN packages through
`install.packages()`. It is safe to run twice. If a package still fails, the
script prints the compiler output; on Debian/Ubuntu the usual missing system
libraries are `libcurl4-openssl-dev libssl-dev libxml2-dev libpng-dev
libfontconfig1-dev libharfbuzz-dev libfribidi-dev`.

Check the installation without starting an analysis:

```bash
Rscript -e 'source("pipeline/run_all.R"); cat("ok\n")'
Rscript pipeline/run_analysis.R --help
```

## 4. Step 2 — run the example

From an R session:

```r
source("pipeline/run_all.R")

run_all(
  counts_file = "input/NIH3T3.xlsx",
  samples     = "input/NIH3T3_metadata.csv",
  species     = "mouse",
  contrast    = c("FBS", "Control"),
  log2fc      = 0.75,     # primary DEG threshold
  min_fpkm    = 5,        # expression floor: higher group-mean FPKM must reach this
  up_cutoff   = 2.5,      # display only: -2.5 to 2.5
  low_cutoff  = 0.15,     # GSEA: minimum absolute true effect; NULL/0 disables it
  top_tfs     = 30        # the setting used in the v1.2 validation
)
```

or from a shell, which runs the same function:

```bash
Rscript pipeline/run_analysis.R \
  --counts input/NIH3T3.xlsx \
  --samples input/NIH3T3_metadata.csv \
  --species mouse \
  --contrast FBS,Control \
  --top-tfs 30 --up-cutoff 2.5 --low-cutoff 0.15
```

The start of the run prints what it understood, and this echo is worth
reading before the analysis gets going:

```
  counts   : input/NIH3T3.xlsx
  samples  : input/NIH3T3_metadata.csv
  species  : mouse (org.Mm.eg.db, mmu, DoRothEA mm)
  expression: max group-mean FPKM >= 5
  contrast : FBS_vs_Control
  thresholds: |log2FC| >= 0.75, adjusted P < 0.05, regulon confidence A/B
  GSEA     : padj < 0.05; |log2FC| >= 0.15; unbounded effects
  display  : log2FC limit +/-2.5
  TFs      : top 30 of the expressed regulons
  sections : 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
  fig sizes: 45 keys from pipeline/R/figure_sizes.R (reloaded for this run)
  output   : results/NIH3T3
```

Output goes to `results/<name of the workbook>/`, so this run writes to
`results/NIH3T3/`. Nothing is overwritten silently: the folder name is fixed,
and `overwrite = FALSE` makes a second run write `NIH3T3_run2` instead.

### Running only part of the workflow

```r
run_all(..., sections = c(1, 2, 3))          # first three sections only
run_all(..., sections = 5, overwrite = TRUE) # redo the TF section after a change
```

Sections that need an upstream result read it from the same output folder, so
re-running section 5 alone reuses the count model of the previous run.

## 5. Step 3 — check that the run is complete

Every section writes its own `status.json`. A quick check:

```bash
for f in results/NIH3T3/section*/status.json; do
  echo "$(dirname "$f"): $(python3 -c "import json,sys;print(json.load(open('$f'))['status'])")"
done
```

Eleven `completed` lines mean the run finished. If a section reports `failed`,
its message names the cause; the rest of the run continues because
`continue_on_error = TRUE`.

The same numbers as the released run, which you can use as a reference:

| Quantity | Value |
| --- | --- |
| genes in the workbook / tested after pre-filter | 55,401 / 14,655 |
| genes above the expression floor (max group-mean FPKM ≥ 5) | 8,369 |
| DEGs at adjusted P < 0.05, \|log2FC\| ≥ 0.75 (primary DEG cut-off) | 1,213 (830 up, 383 down) |
| DEGs at a stricter \|log2FC\| ≥ 1, for reference | 659 (467 up, 192 down) |
| significant GO terms (BP / MF / CC) | 184 / 42 / 63 |
| significant GSEA gene sets (GO BP / MF / CC) | 405 / 61 / 113 |
| significant Hallmark sets | 16 of 46 tested |
| significant KEGG pathways | 27 |
| GSEA ranked identifiers / true log2FC range | 5,412 / -8.459 to 2.766 |
| TF regulons tested / passing the expression floor | 63 / 63 (display limited by `top_tfs = 30`) |
| STRING network | 236 proteins, 1,086 interactions, 10 modules |
| DEG-encoded proteins with GO CC annotation / in selected compartments | 1,100 evaluated / 1,057 mapped across 15 compartments |

The record of *how* the run was made is in `results/NIH3T3/provenance/`:

| File | Content |
| --- | --- |
| `config_used.yml` | every parameter after the defaults were applied — the exact settings |
| `run_manifest.json` | input files, contrast, timestamps, R and package versions |
| `packages.csv` | version of every package that contributed |
| `sessionInfo.txt` | the full R session |
| `figure_sizes.csv` | physical size of every figure drawn |
| `gsea_input_genes.csv` | eligible GSEA input rows with true unshrunk effects and gene padj |
| `gsea_ranked_vector.csv` | exact deduplicated identifier ranking, with no upper clipping |
| `status_gsea_input.json` | GSEA significance policy, lower cutoff and input sizes |

## 6. Step 4 — find your results

```
results/NIH3T3/
├── 00_differential_expression/   DESeq2 model, QC figures and coordinates
├── section01_screening/          DEG tables per threshold, volcano, chromosome landscapes
├── section02_go_ora/             GO over-representation (BP/MF/CC), classes, term–gene networks
├── section03_gsea_go/            GO GSEA, MSigDB Hallmark GSEA, enrichment map
├── section04_kegg_pathview/      KEGG GSEA, lollipop, UpSet, pathview maps
├── section05_tf_regulon/         regulon tables, target enrichment, TF dot plots
├── section06_regulatory_activity/ decoupleR activity, GSVA pathway activity
├── section07_tf_network/         TF–target network, node and edge tables
├── section08_tf_modules/         TF modules, membership and similarity tables
├── section09_geneset_heatmaps/   GO heat maps, TF regulon dot plots
├── section10_ppi_network/        STRING network, hubs, modules, module GO
├── section11_subcellular_localization/ compartment map and tables
├── go_term_target_network/       the multi-functional term–gene network
├── provenance/                   the record described above
└── report/                       figure atlas input, HTML draft, results_report.docx
```

Two conventions make the files easy to find:

- file names carry the section prefix (`S1_`, `TF_`, …) and, where a file
  belongs to a comparison, the contrast name (`…_FBS_vs_Control.csv`);
- every `.csv` is plain, header-first UTF-8 and can be opened in R, Python or a
  spreadsheet without renaming columns.

## 7. Step 5 — answer a question from the files

| Question | Read |
| --- | --- |
| Did the samples behave? | `00_differential_expression/qc_pca_coordinates.csv`, `qc_pca.png`, `qc_sample_correlation.png` |
| Which genes changed? | `section01_screening/S1_DEG_table_lfc_0.75.csv` (`GenesFoldChange` and `GenesFoldChange_raw` retain true unshrunk effects; figures alone are clipped) |
| How robust is the call to the threshold? | `S1_DEG_table_lfc_0.15.csv` … `_2.50.csv`, or `S1_threshold_counts.png` |
| Where are they on the genome? | `S1_chromosome_positions_lfc0.75_FBS_vs_Control.csv` |
| What processes are over-represented? | `section02_go_ora/ENRICHGOBPs_FBS_vs_Control.csv` (simplified) and `ENRICHGOBPalls_…csv` (unfiltered) |
| What is enriched in the independently pre-screened ranking? | `section03_gsea_go/GSEA_GO_BP_FBS_vs_Control.csv`, `GSEA_MSigDB_H_FBS_vs_Control.csv` |
| Which pathways? | `section04_kegg_pathview/KEGG_FBS_vs_Control.csv` |
| Which regulators explain the response? | `section05_tf_regulon/TF_regulon_summary_FBS_vs_Control.csv` — `targets_deg`, `NES`, `gsea_padj`, `coherence`, `tf_log2FC` |
| Is the regulator active, independent of the DEG list? | `section06_regulatory_activity/TF_activity_diff_ulm.csv` |
| Who regulates whom? | `section07_tf_network/TF_network_edges_FBS_vs_Control.csv`, `TF_network_nodes_…csv` |
| Which TFs share targets? | `section08_tf_modules/TF_cluster_modules_FBS_vs_Control.csv` |
| Which proteins interact? | `section10_ppi_network/PPI_edges_FBS_vs_Control.csv`, `PPI_nodes_…csv` (degree, betweenness, module), `PPI_modules_…csv` |
| Where do the proteins sit? | `section11_subcellular_localization/subcellular_localization_FBS_vs_Control.csv` |

A caution that the workflow is built around: a regulon with many DE targets and
a large NES but a low `coherence` means the target set overlaps the response
without the TF following it. TP53 is the shipped example — enriched targets,
falling inferred activity.

## 8. Step 6 — change the analysis

Most day-to-day changes are single arguments:

```r
# stricter primary DEG calls (GSEA retains its independent lower bound)
run_all(..., log2fc = 1, alpha = 0.01)

# change gene display scales only; calculation tables retain true effects
run_all(..., up_cutoff = 3)

# remove the GSEA effect-size lower bound, keeping gene padj < alpha
run_all(..., low_cutoff = NULL) # low_cutoff = 0 is equivalent

# a stricter expression floor: only genes reaching 10 FPKM in one group
run_all(..., min_fpkm = 10)

# more transcription factors, all DoRothEA tiers
run_all(..., top_tfs = 60, confidence = c("A", "B", "C"))

# no KEGG downloads, no STRING query
run_all(..., pathview_download = FALSE, run_ppi = FALSE)

# human data
run_all(counts_file = "input/human.xlsx", species = "human")
```

Every argument and its default is listed in the README table. If you change one
of the analysis thresholds, delete nothing: re-run with `overwrite = TRUE` (or
into a new folder) so the provenance record matches the figures.

Two settings deserve a warning because they change meaning rather than
strictness:

- `gene_length_column` — the workbook column used for FPKM. If your workbook
  has no length column, set it to `NULL`; the pipeline then uses a GTF
  (`gene_annotation`), a `*.gtf` beside the workbook, or the species
  annotation package.
- `expression_matrix_file` — supply a pre-computed FPKM/TPM matrix when you
  have one; the heat-map sections use it directly instead of computing FPKM.

## 9. Using your own data

Before the first run with a new experiment, check four things:

1. **One row per gene, one column per library**, plus an identifier column and,
   if you have them, symbol and description columns. Excel and CSV both work.
2. **Raw integer counts** for the count columns — the model is a negative
   binomial on counts; FPKM or TPM in those columns is not a substitute.
3. **A sample table** with `sample_id`, `condition` and (optionally)
   `replicate`, matching the count columns exactly.
4. **An expression floor you are happy with.** By default a gene must reach
   5 FPKM (group mean) in at least one arm to enter the analysis; set
   `min_fpkm` to `NULL` to analyse every tested gene.
5. **Two groups only.** More than two levels are accepted but the workflow
   analyses one contrast at a time (`contrast = c("A", "B")`).

If the workbook carries an exon-union length column, name it in
`gene_length_column`; otherwise let the pipeline take lengths from the
annotation. If you want to see the length table and the computed FPKM, they are
written to `00_differential_expression/gene_lengths_used.csv` and
`FPKM_computed.csv`.

## 10. When something goes wrong

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| `Metadata sample_id values do not exactly match count-matrix columns` | typo, or the workbook has an extra column (e.g. a summary column) | make the six strings identical; pass `count_columns` to name the sample columns explicitly |
| `Required package(s) missing` | installation incomplete | `Rscript scripts/install_packages.R`, then re-run |
| Section 4 empty or slow | KEGG download blocked | the run falls back to the cache; set `pathview_download = FALSE` to skip the maps and keep the statistics |
| Section 10 `STRING request failed` | no network | the request is retried up to `ppi_timeout`; run again later or with `run_ppi = FALSE` |
| Section 9 empty | no expression matrix and no TF resource | check `section09_geneset_heatmaps/status.json`; FPKM falls back to normalised counts, and the TF part needs `section05` |
| `results_report.docx` missing | the fixed-template report failed after the analysis | see `report/status_summary.json`; the analysis outputs are unaffected |
| Figures look different from the paper | a parameter differs | compare `provenance/config_used.yml` with the settings in the README table |

Section failures never take the whole run down: each section records its own
status, and `continue_on_error = FALSE` is available when you would rather stop
at the first failure.

## 11. Rebuilding the paper

The figures and tables of the accompanying paper are rebuilt from the run, not
hand-edited:

```bash
python3 work/make_paper_tables.py --results results/NIH3T3   # Tables 1–7
python3 work/build_word_figures.py                           # composite figures
Rscript scripts/build_figure_atlas.R results/NIH3T3          # figure atlas PDF
```

`work/README.md` describes the full chain to the Word and PDF manuscripts.

## 12. Citing and asking for help

Cite the software with the record in `CITATION.cff` (and the paper, when it is
published). Questions and bug reports: open an issue at
<https://github.com/ChangyeSun/TransAtlas-R/issues>. When reporting a failed
section, attach its `status.json` and the matching console output — between
them they usually contain the answer.

---

Author: Changye Sun, Henan Key Laboratory of Tissue Regeneration, Henan Medical
University, Xinxiang, Henan 453003, China · sunchangye-boy@163.com.
Code organisation and documentation were prepared with the assistance of
ClawsGO (ClawsGO Science Agent).
