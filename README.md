# TransAtlas-R — generic two-group bulk RNA-seq workflow

**Version 1.2.0** · [MIT licence](LICENSE) · [changelog](CHANGELOG.md) · [citation](CITATION.cff)

**Author:** Changye Sun — Henan Key Laboratory of Tissue Regeneration, Henan Medical University, Xinxiang, Henan 453003, China · sunchangye-boy@163.com

Changye Sun conceived the workflow, implemented the analysis sections, designed and ran the case study, and wrote the manuscript. ClawsGO assisted with code organisation and documentation.

**TransAtlas-R** is the project name: one entry point, two input files, eleven paper-ordered analysis sections, and automated Word/HTML/atlas deliverables.

This project contains a configurable, publication-oriented **two-group bulk RNA-seq** workflow. It is designed for a general experimental group versus control group and is independent of any particular organism, tissue, cell line, treatment, disease model, or molecular pathway.

The workflow covers differential expression, functional enrichment, transcription-factor regulation, protein interaction and subcellular localization. The released example is the six-sample NIH 3T3 serum-restimulation experiment of the accompanying paper.

Run the commands below from the project root. Standalone plotting and atlas scripts also use this working-directory convention.

## Installation

**Requirements.** R 4.2 or newer (developed on R 4.6.1 / Bioconductor 3.23), Pandoc 2.11 or newer for the HTML and Word report, and a C/C++ toolchain for packages built from source. No GPU is needed.

```bash
git clone https://github.com/ChangyeSun/TransAtlas-R.git
cd TransAtlas-R

Rscript scripts/install_packages.R                    # mouse and shared packages
Rscript scripts/install_packages.R --species human    # add the human annotation set
```

The script installs exactly the packages `run_all()` checks at start-up, so the two lists cannot drift apart: Bioconductor packages (DESeq2, clusterProfiler, fgsea, GSVA, decoupleR, dorothea, pathview, the OrgDb/EnsDb annotation packages) go through `BiocManager`, CRAN packages through `install.packages()`. It is safe to repeat — installed packages are skipped — and it prints the build log of anything that still fails.

System libraries that a source installation usually needs:

```bash
# Debian / Ubuntu
sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libpng-dev \
                 libfontconfig1-dev libharfbuzz-dev libfribidi-dev pandoc
# macOS
brew install pandoc libxml2 cairo pkg-config
```

Nothing else has to be downloaded in advance. The DoRothEA regulons for mouse and human ship in `data/dorothea/` (the mouse collection also serves rat); the KEGG pathway definitions and STRING responses are cached under `cache/` on first use, so the first sections 4 and 10 need network access and later runs do not.

## Sample run

A run is one call to `run_all()`. Everything is an argument, nothing is read from a configuration file, and anything omitted keeps the built-in default.

```r
source("pipeline/run_all.R")

run_all(
  counts_file = "input/NIH3T3.xlsx",   # omit -> file-selection dialog
  species     = "mouse",               # mouse | human | rat
  log2fc      = 0.75,                  # the one |log2FC| cut-off of the workflow
  min_fpkm    = 5,                     # expression floor: max group-mean FPKM
  top_tfs     = 30                     # transcription factors analysed
)
# -> results/NIH3T3/
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `counts_file` | `NULL` → file dialog | raw-count matrix (.xlsx/.xls/.csv/.txt); the output folder takes its name |
| `samples` | `NULL` | metadata CSV/Excel/TSV or data.frame (`sample_id`, `condition`, `replicate`); `NULL` finds `<data>_metadata.csv` / `sample_metadata.csv` next to the data file, then in `config/`, and accepts it only when the sample ids match the count columns |
| `species` | `"mouse"` | `"mouse"`, `"human"` or `"rat"`; sets the OrgDb, KEGG code, identifier types, MSigDB database and DoRothEA collection (hs / mm) in one word |
| `group_column` | `"condition"` | metadata column holding the two groups |
| `contrast` | `NULL` | `c(numerator, denominator)`, e.g. `c("Treatment", "Control")`; `NULL` uses the metadata levels (second alphabetically over the first) |
| `log2fc` | `0.75` | primary absolute log2FC threshold for DEG calls and DEG gene sets; GSEA uses the independent `low_cutoff` |
| `min_fpkm` | `5` | expression floor shared by every section: the higher of a gene's two group-mean FPKM values must reach this value, otherwise the gene is dropped before screening, enrichment and the regulatory sections; `NULL` switches it off |
| `alpha` | `0.05` | adjusted-P cut-off |
| `min_count` | `10` | minimum raw count in `min_samples` samples for a gene to be tested |
| `min_samples` | `2` | number of samples that must reach `min_count` |
| `gene_id_column`, `gene_symbol_column`, `annotation_column` | `"gene_id"`, `"GeneSymbol"`, `"Description"` | annotation columns of the count workbook |
| `gene_length_column` | `"Exonic.gene.sizes"` | workbook column holding the exon-union length of each gene; it is the primary source for the computed FPKM. `NULL` or `"none"` disables it |
| `counts_sheet` | `NULL` | Excel sheet; `NULL` = first sheet |
| `count_columns` | `NULL` | sample columns of the count matrix; `NULL` infers them from the metadata |
| `expression_matrix_file` | `NULL` | optional pre-computed FPKM/TPM matrix; `NULL` first looks for `<counts name>_fpkm.csv` beside the count workbook, then computes FPKM from the raw counts and the exon-union lengths (the workbook column first, the annotation when the column is absent), and records the length source and identifier coverage; missing FPKM stops a run with an active expression floor |
| `expression_matrix_unit` | `"FPKM"` | unit of a supplied matrix: `"FPKM"`, `"TPM"` or `"normalized_counts"`; TPM/normalized counts require `min_fpkm = NULL` |
| `gene_annotation` | `NULL` | GTF/GFF used for the exon-union lengths when the workbook has no length column; `NULL` searches `options(transatlas.gene_annotation=)`, then `TRANSATLAS_GENE_ANNOTATION`, then a `*.gtf` beside the count workbook, and finally the species annotation package (`EnsDb.Mmusculus.v79`, `EnsDb.Hsapiens.v86`, `TxDb.Rnorvegicus.UCSC.rn6.refGene`) |
| `top_tfs` | `30` | number of transcription factors analysed and drawn (dot plot and network) |
| `confidence` | `c("A","B")` | DoRothEA confidence tiers kept when building regulons |
| `tf_min_targets`, `tf_ora_min_overlap`, `tf_coherence_min_targets` | `5`, `3`, `5` | regulon size floor, minimum DEGs in a regulon for ORA, minimum targets for the coherence score |
| `tf_file` | `NULL` | optional custom TF–target table; `NULL` uses the DoRothEA data for the chosen species |
| `tf_dorothea_dir` | `NULL` → `data/dorothea` | folder with the bundled DoRothEA `.rda` files |
| `max_genes` | `30` | upper bound on genes drawn in any heat map or dot plot |
| `top_n_labels` | `20` | labels drawn on the volcano / effect panels |
| `volcano_label_expression_n`, `volcano_label_effect_n` | `5`, `5` | expression-ranked and effect-ranked labels per direction, within `top_n_labels` |
| `up_cutoff` | `2.5` | display-only symmetric log2FC limit (−2.5 to 2.5); `NULL` uses an automatic display range. All calculation tables retain true values |
| `low_cutoff` | `0.15` | independent GSEA minimum absolute log2FC; `NULL` or `0` disables the lower bound, retaining the gene `padj < alpha` filter |
| `log2fc_cap` | `NULL` | deprecated display-only alias of `up_cutoff` |
| `volcano_lfc_thresholds` | `c(0.15, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5)` | threshold ladder for screening tables and figures; independent of GSEA input |
| `expression_log10_threshold` | `NULL` | legacy count-based expression rule (`log10(baseMean) <= 2`); passed `2` it hides those genes from the up/down classes, which conflicts with the single FPKM floor, so it is off by default |
| `go_pvalue_cutoff`, `go_qvalue_cutoff` | `0.05`, `0.05` | GO over-representation cut-offs |
| `go_bp_pvalue_cutoff`, `go_bp_qvalue_cutoff` | `0.01`, `0.01` | stricter cut-offs applied to GO BP |
| `go_simplify_cutoff` | `0.5` | semantic-similarity cutoff for GO simplification |
| `gsea_min_size`, `gsea_max_size`, `gsea_pvalue_cutoff` | `15`, `500`, `0.05` | GSEA size window and reporting cutoff |
| `gsea_top_sets` | `9` | sets drawn on a composite GSEA panel |
| `run_msigdb`, `msigdb_category` | `TRUE`, `"H"` | MSigDB collection GSEA |
| `pathview_download`, `pathview_max_pathways` | `TRUE`, `12` | draw pathview KEGG pathway maps; `pathview_download = FALSE` skips the KGML/native-map downloads and the maps, while the rest of section 4 still runs. `run_pathview` remains as a deprecated alias that takes precedence when supplied |
| `network_layouts` | `c("fr","kk","stress","circle")` | layouts rendered for the TF–target network |
| `network_max_labels` | `30` | target labels in the network |
| `cluster_min_jaccard`, `cluster_max` | `0.1`, `8` | section 8 links two TFs when their DE target sets share at least this Jaccard index, and cuts at most this many modules |
| `activity_methods` | `c("ulm","aucell","wmean")` | decoupleR methods for regulon activity |
| `activity_min_size`, `activity_alpha` | `5`, `NULL` → `alpha` | regulon size floor and significance level for TF/GSVA activity; an explicit value overrides `alpha` |
| `run_gsva`, `gsva_collection`, `gsva_min_size`, `gsva_max_size` | `TRUE`, `"H"`, `10`, `500` | GSVA pathway activity |
| `run_ppi`, `ppi_score_threshold`, `ppi_max_nodes`, `ppi_label_n`, `ppi_max_modules` | `TRUE`, `0.4`, `300`, `4`, `10` | STRING protein–protein interaction network of the DEG-encoded proteins; confidence threshold, maximum proteins submitted, maximum hub labels per module, and the module ceiling (smallest communities are merged into their best-connected neighbour, and every module is annotated with GO terms) |
| `ppi_module_annotation_alpha` | `0.2` | exploratory module GO annotation cutoff; output marks statistical significance separately using `alpha` |
| `ppi_cache_dir`, `ppi_timeout` | `NULL` → `cache/string`, `120` | where STRING responses are cached and the request timeout |
| `qc` | `TRUE` | draw the quality-control figures |
| `qc_plots` | `NULL` → all | subset such as `c("pca","correlation","ma_plot","dispersion","library_size","count_distribution","independent_filter","sample_distance","rlog_vst")` |
| `sections` | `NULL` | e.g. `c(1, 3, 6)`; `NULL` runs the default set |
| `report` | `TRUE` | render the HTML report at the end of the run |
| `report_title`, `report_author`, `report_language`, `generate_pdf` | see `run_all()` | report metadata |
| `kegg_cache_dir`, `kegg_download_timeout` | `NULL`, `240` | where KEGG tables are cached and the first-download timeout |
| `check_packages`, `install_missing_packages` | `TRUE`, `TRUE` | check every package before the run, install what is missing, re-check, then analyse; a missing required package stops the run with the exact list |
| `continue_on_error` | `TRUE` | keep going when an optional module fails |
| `output_dir` | `NULL` | defaults to `results/<name of the data file>/` |
| `overwrite` | `TRUE` | `FALSE` writes to `<folder>_run2`, `_run3`, … instead of reusing the folder |
| `seed` | `20261013` | random seed for the reproducible steps |
| `config_file` | `NULL` | optional YAML overlay; not needed for a normal run |

Only two things have to be supplied: the count matrix and the sample metadata. TF–target regulons are read from the DoRothEA data distributed with the `dorothea` package (`dorothea_hs` for human, `dorothea_mm` for mouse and rat); when the package is not installed, the identical `.rda` objects in `data/dorothea/` are used, so the workflow runs offline.

The same settings from the command line:

```bash
Rscript pipeline/run_analysis.R --counts input/NIH3T3.xlsx \
        --species mouse --log2fc 0.75 --min-fpkm 5 --top-tfs 30
# re-run selected sections only
Rscript pipeline/run_analysis.R --counts data.xlsx --species human --sections 1,3,6
# explicit contrast and wider regulons
Rscript pipeline/run_analysis.R --counts data.xlsx --contrast Treatment,Control \
        --confidence A,B,C --max-genes 40
```

The CLI mirrors every `run_all()` argument: input and column names, `--log2fc`, `--alpha`, `--min-count`, `--min-samples`, `--group-column`, `--contrast A,B`, `--top-tfs`, `--confidence A,B`, `--min-fpkm`, `--max-genes`, all enrichment cut-offs, `--network-layouts`, `--activity-methods`, `--qc-plots`, `--sections`, `--no-report`, `--no-pathview`, `--no-pathview-download`, `--no-gsva`, `--no-msigdb`, `--no-qc`, `--out` and `--keep-existing`. The header of the script lists every flag.

Results are written to `results/<data file name>/`. A DESeq2 fit is cached under `<output_dir>/.cache/` and reused only when the selected counts, annotations, complete metadata, design, contrasts, filtering, alpha, shrinkage, seed and relevant package versions match, so individual sections can be re-run quickly. Sections are independently runnable: the shared screened-gene context is rebuilt when a run starts after section 1, and sections that consume another section's table reload it from the output directory (section 7 reads section 5's table, section 9 reads the GO result of section 2).

The run writes `run_status.json` with `completed`, `completed_with_errors` or `failed`. A screening failure stops the run. Other section failures preserve the last valid context when `continue_on_error = TRUE`; the CLI exits nonzero after recording failed steps. Optional operation errors are retained in `optional_errors.log` within the relevant section. Each run also writes `warnings.csv` (stage, call, class, message, UTC time and run ID), `warning_summary.csv` and `status_warnings.json`; warnings remain visible to R. Stage folders retain their own warning records. Report warnings are included in the same central record. These files describe the latest invocation, and are finalized even when the run fails. Activity tables retain profiles with undefined Welch tests using `test_status`, `test_reason` and NA P values; their descriptive statistics remain available and BH correction includes only valid tests. GSEA and GSVA write `*_source.json` with the actual MSigDB version/database, collection, target species, ortholog policy and a checksum of the gene-set membership. Missing significant GO terms are a valid empty result, and non-significant fallback terms are never reloaded as significant hits.

## Sections

Each file in `pipeline/sections/` implements one analysis section.

| # | Section | Content |
| --- | --- | --- |
| 1 | `s01_screening.R` | Multi-threshold DEG tables and `sData.xlsx`, five-way gene classification, volcano and expression–effect panels, threshold ladder, chromosome landscapes at \|log2FC\| ≥ 0.5, 0.75 and 1.0, global and top-DEG heat maps |
| 2 | `s02_go_ora.R` | GO over-representation per ontology, semantic simplification, dot plots, rich-factor plot, functional-class map and term–gene networks |
| 3 | `s03_gsea_go.R` | Ranked GSEA for GO BP/MF/CC and MSigDB Hallmark, composite panels, enrichment map and one-page-per-set PDFs |
| 4 | `s04_kegg_pathview.R` | KEGG enrichment from the local pathway cache, lollipop and UpSet plots, multi-pathway GSEA panel and optional pathview maps |
| 5 | `s05_tf_regulon.R` | TF–target edges, regulon over-representation, regulon GSEA, coherence score, TF dot plots and chromosome landscapes of the strongest regulons |
| 6 | `s06_regulatory_advanced.R` | Regulon activity (decoupleR ULM/AUCell/WMean) and GSVA pathway activity |
| 7 | `s07_tf_network.R` | TF–target network in force-directed, Kamada–Kawai, stress and circular layouts |
| 8 | `s08_tf_clusters.R` | TF modules from the Jaccard overlap of their DE target sets, plus the module network |
| 9 | `s09_geneset_heatmaps.R` | GO term heat maps and TF-regulon dot plots |
| 10 | `s10_ppi.R` | STRING protein–protein interaction network, hub ranking and module GO annotation |
| 11 | `s11_subcellular.R` | GO cellular-component mapping and compartment enrichment |

Sections are numbered in the order a paper reports them: differential expression → GO over-representation → GSEA → KEGG → TF regulation (regulon, activity, network, modules) → gene-set display → protein interaction → subcellular localization.

Shared infrastructure lives in `pipeline/R/`: `helpers.R` (configuration, input reading, plotting defaults), `context.R` (the per-run context passed between sections, including the FPKM matrix and the TF resource), `deseq2_module.R`, `qc_module.R`, `provenance.R`, and `report_module.R`.

## Generic input contract

### Count matrix

Required columns:

- one stable gene identifier column, e.g. `gene_id`;
- one optional gene-symbol column, e.g. `GeneSymbol`;
- one integer raw-count column per sample;
- optional annotation columns such as a gene description, and an optional `Exonic.gene.sizes` column with the exon-union gene length used for FPKM.

The pipeline accepts Excel, CSV, and tab-delimited text. It rejects negative, non-integer, missing, or non-numeric count values, and sums duplicate stable identifiers while keeping the first annotation.

### Expression matrix for the heat-map sections

TransAtlas takes that matrix as `expression_matrix_file`; when it is left `NULL` the workflow looks for `<counts name>_fpkm.csv` next to the count workbook, and if there is none it computes FPKM itself from the raw counts and the exon-union gene lengths. The lengths come first from the workbook column named by `gene_length_column` (`Exonic.gene.sizes` by default, as in `input/NIH3T3.xlsx`); when the workbook carries no such column they come from a GENCODE or Ensembl GTF — pass it as `gene_annotation = "gencode.vM25.annotation.gtf"` (or `--gene-annotation` on the command line), or set it with `options(transatlas.gene_annotation = ...)` or the `TRANSATLAS_GENE_ANNOTATION` environment variable, or simply place a `.gtf` next to the count workbook — and otherwise from the species annotation package (`EnsDb.Mmusculus.v79`, `EnsDb.Hsapiens.v86`, `TxDb.Rnorvegicus.UCSC.rn6.refGene`). A missing workbook length column therefore still permits FPKM calculation through the annotation database. Individual genes without matching lengths are omitted and counted in `status_fpkm.json`. If every length source fails, a run with an active `min_fpkm` stops with an explicit diagnostic. Setting `min_fpkm = NULL` explicitly allows display fallback to DESeq2 normalized counts. The length table used and the computed matrix are written to `00_differential_expression/gene_lengths_used.csv` and `FPKM_computed.csv`; annotation tables are cached under `cache/gene_lengths/`. FPKM values define the shared expression floor and display matrix; DESeq2 always models raw counts.

### Sample metadata

```text
sample_id,condition,replicate
Control-1,Control,1
Control-2,Control,2
Treatment-1,Treatment,1
Treatment-2,Treatment,2
```

CSV, Excel (`.xlsx`/`.xls`), TSV and in-memory data.frames use the same metadata reader. Automatic discovery accepts exact sample matches and skips subset or duplicated sample lists. The group column is configurable (`group_column`). The workflow requires exactly two levels and checks the numerator/denominator labels before modelling. Biological replicates must be independent experimental units.

## Defaults and overrides

There is no configuration file. Every default lives in `default_config()` in `pipeline/R/helpers.R`, and every value can be overridden by the corresponding `run_all()` argument listed under [Sample run](#sample-run). A YAML overlay is still accepted through `run_all(config_file = "path/to/settings.yml")` or `run_analysis.R --config ...` for users who prefer one, but nothing in a normal run reads it.

### Species selection

`species` takes one word and replaces the whole annotation block:

| `species` | OrgDb | KEGG | Ensembl prefix | MSigDB database | DoRothEA |
| --- | --- | --- | --- | --- | --- |
| `"mouse"` | `org.Mm.eg.db` | `mmu` | `ENSMUSG` | `MM` | `mm` |
| `"human"` | `org.Hs.eg.db` | `hsa` | `ENSG` | `HS` | `hs` |
| `"rat"` | `org.Rn.eg.db` | `rno` | `ENSRNOG` | `RN` | `mm` |

Install and verify the matching OrgDb package. The preset also sets `kegg_key_type: ncbi-geneid`, which clusterProfiler's KEGG interface expects for Entrez identifiers, and `pathview_id_type: ENTREZ`, which is pathview's own vocabulary. The module automatically falls back to the ortholog-mapped human MSigDB collections when the requested collection is unavailable (`H` requests human Hallmark, including ortholog mapping for mouse/rat; `MH` explicitly requests mouse-native Hallmark where available). Collection availability is checked before retrieval, and the exact gene-set source is recorded for both GSEA and GSVA.

### Transcription-factor resource

Transcription-factor regulons come from the **DoRothEA** collection distributed with the `dorothea` package: `dorothea_hs` for human and `dorothea_mm` for mouse. DoRothEA publishes no rat-specific set and rat symbols follow mouse spelling, so `species = "rat"` uses `dorothea_mm`. `confidence` (default `c("A", "B")`) selects the confidence tiers kept when regulons are built.

When the `dorothea` package is not installed, the identical `.rda` objects bundled in `data/dorothea/` are loaded instead, so a run needs no TF file and no network access. `tf_file` remains as an override: give any table with `TF`, `target` and optionally `confidence` / `mor` columns (for example CollecTRI or TRRUST) and it replaces the DoRothEA data for that run. When no resource can be read, the run warns with the expected path and TF-dependent sections record a skipped status.

TF–target networks are descriptive hypotheses. They do not establish direct binding or causal regulatory activity; that requires ChIP-seq, CUT&RUN, reporter assays, or perturbation experiments.

### Choices that affect the results

- **One expression floor.** A gene enters the analysis when the higher of its two group-mean FPKM values reaches `min_fpkm` (default 5); genes below it are dropped before screening, enrichment and the regulatory sections, so every section works on the same expressed universe. The floor is written to `00_differential_expression/expression_floor.csv` (`max_group_mean_fpkm`, `kept`); when FPKM cannot be built, an active floor stops the run and `status_expression_floor.json` records why; `min_fpkm = NULL` disables it explicitly. The transcription-factor figures use the same floor.
- **True effects are retained for every calculation and result table.** `up_cutoff = 2.5` clips only gene display coordinates, stem heights and colour scales. `GenesFoldChange` and `GenesFoldChange_raw` both contain the unbounded unshrunk estimate. `low_cutoff = 0.15` independently defines the GSEA lower bound, together with gene `padj < alpha`. `low_cutoff = NULL` or `0` removes the effect-size lower bound. The DEG threshold remains `log2fc = 0.75`. ORA uses the eligible tested universe. GSEA input genes and the exact ranked vector are exported in `provenance/`.
- **GSEA panels are non-redundant and two-sided.** Sets with more than 50 genes are ranked by \|NES\|; a set sharing more than 40% of the smaller set's genes with one already chosen is skipped; up/down slots follow the positive/negative NES counts. Result tables keep every set.
- **KEGG data are cached.** The pathway table is fetched once into `cache/kegg/` and reused. If it cannot be fetched, section 4 falls back to the MSigDB KEGG collection and records the source in `KEGG_source_<contrast>.txt`; pathview maps are skipped because they need real KEGG ids. `pathview_download = FALSE` skips only the map downloads.
- **Runs are serial.** The pipeline registers `BiocParallel::SerialParam()` and sets `mc.cores = 1` to keep memory bounded on shared hosts.

## Figure style

The house style is set in one place, `pipeline/R/helpers.R`, and every figure follows it.

- **Heat maps** share one diverging ramp, `heatmap_palette()`: purple for the low end through a light middle to orange for the high end (ColorBrewer PuOr). It replaced the green-white-blue ramp, whose two ends both read as "blue-green" once pheatmap stretched the row z-scores over the full cell range. Applied to the screening heat maps, the gene-set heat maps and the TF-activity heat maps.
- **Functional-class map.** The GO terms clustered into ten classes are drawn as a coloured network in a saturated qualitative palette with darker edges; the class names live in the legend only, so no label sits on top of the class it names. Legend text and keys are 50% larger than in earlier versions (19.5/21 pt, 27 pt keys), and the ten class symbols are two steps larger than the default key glyph.
- **Edges** in the TF-target network carry the mode of regulation in saturated red (`activating`) and blue (`repressing`) at half opacity, which stays legible where several hundred edges overlap.
- Base font sizes are 13 pt in the heat maps, 12 pt in figure legends and 15 pt in titles.
- Size legends use `size_breaks()`: breaks are clipped to the data range and always include its minimum and maximum, so the smallest and largest keys match the smallest and largest dots instead of a `pretty()` endpoint that lies beyond the panel.

## Figure sizes and fast redraws

Every saved figure takes its width and height from a single registry, `pipeline/R/figure_sizes.R` (`FIGURE_SIZES`), in inches. No section hard-codes a size any more. An entry is either `c(width, height)` or a function of the rows shown, for the panels whose height follows their content (`s2_go_dotplot`, `s6_tf_dotplot`, `s6_tf_expression`).

The registry is re-read from disk at the start of **every** `run_all()` call (and of every direct `run_pipeline()` call), so an edited number applies to the next run in the same session - no re-sourcing of `run_all.R` and no restart. Each run then writes the sizes it actually used to `provenance/figure_sizes.csv` (key, width, height, file), one row per figure, so the dimensions behind any published panel can be checked after the fact.

To change a figure size:

1. Edit the number in `pipeline/R/figure_sizes.R`.
2. Re-save the affected figures **without re-running the analysis**:

```bash
# every cached ggplot figure
Rscript scripts/redraw_figures.R

# one figure by registry key, or several
Rscript scripts/redraw_figures.R s3_go_composite
Rscript scripts/redraw_figures.R s1_volcano s8_tf_network

# what is available, with the current size
Rscript scripts/redraw_figures.R --list
```

`save_plot()` stores each ggplot object in `cache/figures/<figure file name>.rds`, so the redraw reads the plot that was already computed and only writes it again at the new size. This covers every ggplot figure, including the ggraph networks and the patchwork composite panels.

Four kinds of figure are drawn directly by their device (pheatmap, base graphics, UpSetR, `gseaplot2` panels): the QC heat maps, the two screening heat maps, the gene-set heat maps, the TF-activity heat map, the dispersion/MA plots, the GSEA multi-pathway panel and the UpSet plot. Their sizes also come from the registry, but there is no plot object to replay, so after the edit rerun just the owning section — this costs seconds to a couple of minutes, not a full `run_all()`:

```bash
Rscript pipeline/run_analysis.R --counts input/NIH3T3.xlsx --species mouse \
        --sections 3 --no-report      # GSEA panel
Rscript pipeline/run_analysis.R --counts input/NIH3T3.xlsx --species mouse \
        --sections 1 --no-report      # screening heat maps / volcano PDF
```

The section numbers are 1–11; `--no-report` skips the manuscript render. The same thing from R:

```r
source("pipeline/run_all.R")
run_all(counts_file = "input/NIH3T3.xlsx", species = "mouse",
        sections = c(3, 6, 7), report = FALSE)
```

After a redraw, rebuild the PDF atlas if you need it to match:

```bash
Rscript scripts/build_figure_atlas.R
```

The redraw only rewrites PNGs that already exist; a figure whose section has never run has no cache and must be produced once by running that section.

The cache lives in `cache/figures/` (about 0.5 GB for the demo data, mostly the 18-inch network panels). It is disposable: delete it at any time and the next section run recreates the entries it needs. Set `CLAWSGO_FIGURE_CACHE` to put it elsewhere.

## Standalone GO term-target network

`scripts/go_term_gene_network.R` draws a bipartite term–gene network for any list of GO terms, independently of the pipeline run. Each requested term becomes a square node, each annotated gene a circle coloured by its log2 fold change, and the layout is Kamada–Kawai (`--layout kk`; `fr`, `stress` and `circle` are also accepted). A term's targets are its own and its descendants' annotations, which is the gene set an enrichment test scores.

```bash
Rscript scripts/go_term_gene_network.R \
    --terms "GO:0006281,mitotic cell cycle process,extracellular matrix organization" \
    --results results/NIH3T3 --species mouse \
    --out results/NIH3T3/go_term_target_network
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--terms` | — | comma- or semicolon-separated GO ids or exact term names |
| `--term-file` | — | text file with one id or name per line, added to `--terms` |
| `--de` | unshrunk `DESeq2_*.csv` under `--results`, with shrinkage table as fallback | differential-expression table to read |
| `--up-cutoff` | `2.5` | display-only symmetric log2FC colour limit; exported effects remain unbounded |
| `--species` | `mouse` | `mouse`, `human` or `rat` (selects the OrgDb) |
| `--out` | `<results>/go_term_target_network` | output directory |
| `--max-genes` | `80` | genes drawn per term, strongest \|log2FC\| first (0 = all) |
| `--label-top` | `40` | gene labels drawn, strongest \|log2FC\| first |
| `--layout` | `kk` | `kk`, `fr`, `stress` or `circle` |
| `--width`, `--height` | `14`, `10` | canvas in inches |

It writes `GO_term_target_network.png` and `.pdf` plus `GO_term_target_edges.csv`, `GO_term_target_nodes.csv` and `GO_term_target_summary.csv` (genes in term, genes drawn, strongest target). `--help` prints the same list.

## Fixed-template result report

Every run also writes `report/results_report.docx`, an English Word report whose thirteen sections never change order or field names. The template lives in `pipeline/R/summary_report.R` as `SUMMARY_TEMPLATE`; after the analysis finishes, `render_summary_report()` reads the CSV tables the sections already wrote and fills each section with that run's numbers.

| Section | Filled with |
| --- | --- |
| 1 Data and design | data file, species, group sizes, contrast, thresholds, genes tested, DEG counts |
| 2 Differential expression | DEG counts at every \|log2FC\| cut-off, strongest up/down genes |
| 3 GO functions | significant terms per ontology, ten leading terms per ontology, the ten functional classes |
| 4 GSEA | significant sets and up/down split, strongest activated and repressed GO / MSigDB sets |
| 5 KEGG | significant pathways with direction/NES, top pathways |
| 6 TF regulons | tested/significant/expressed TFs, leading TFs by target count, key targets |
| 7 Activity | ULM/AUCell/WMean differential activity, GSVA pathway activity |
| 8 TF network | TF, target and edge counts, activating/repressing edges, TF list |
| 9 TF modules | module members, target union, representative TF and GO BP |
| 10 Gene-set display | the 10 GO terms as heat maps and the 10 TF regulons as dot plots |
| 11 PPI | STRING nodes/edges, hubs, every module function |
| 12 Localization | compartment counts and enrichment |
| 13 Conclusions | eleven data-driven statements summarising the sections above |

The same run writes `report/summary_report_values.csv` (section, table, item, value) so the numbers can be reused downstream. The appendix embeds eight key figures at their native aspect ratio. Set `summary_report = FALSE` (`--no-summary-report`) to skip it, or regenerate it without re-running the analysis:

```bash
Rscript scripts/build_summary_report.R \
    --results results/NIH3T3 \
    --out results/NIH3T3/report
```

`--config` points at another `config_used.yml`; `--file` renames the output.

### If `results_report.docx` is missing

Report generation is now a monitored step. A failed run leaves three traces in `report/` instead of ending silently:

- `summary_report_docx_error.txt` — the full R error, platform, library paths and traceback;
- `status_summary.json` — `"status": "failed"` with the same message;
- `summary_report_docx_warnings.txt` — non-fatal warnings collected while writing.

The console also prints `summary_report_docx failed: ...` and the path to `summary_report_docx_error.txt`. The most common cause on macOS is that `officer` or `flextable` is absent or cannot be loaded: `flextable` depends on `systemfonts` and `textshaping`, which need a C++ toolchain. Install Xcode command-line tools once, then the two writers:

```bash
xcode-select --install
Rscript -e 'BiocManager::install(c("officer", "flextable"))'
```

`run_all()` checks both packages at start-up and prints `Word report writers: officer and flextable available` (or `UNAVAILABLE (...)`) before any analysis begins. If the writers cannot be installed, use `summary_report = FALSE` to run the analysis without the DOCX; the CSV tables and figures are still written, and `scripts/build_summary_report.R` can build the report later on a machine where the writers are available.

## Statistical safeguards

- DESeq2 receives raw integer counts only; normalized counts and VST values are used for visualization and correlation-based modules only.
- Significance is based on adjusted P values; effect-size thresholds are reported separately.
- Ranked GSEA uses the independent `low_cutoff` (default 0.15), gene `padj < alpha` and unbounded effects. Changing `up_cutoff` or the screening ladder does not change the ranking.
- Over-representation analysis and TF regulon enrichment both use an explicit tested-gene universe.
- Regulon claims carry an over-representation test, a ranked-list GSEA statistic, and a coherence score, not just a target count.
- Two-group validation and contrast-level checks occur before modelling.
- TF modules are a hypothesis about co-regulated output, not evidence of co-binding: a shared target can be regulated by both factors independently.
- Subcellular compartment counts are GO annotations, not a partition of the DEG list; a protein can map to several compartments, and the expected count comes from all tested genes with the same annotations.
- Optional modules stop transparently when their species resources are unavailable; every section writes a `status.json`.
- Dependencies are checked before any input is read (`check_packages`); the check covers the analysis sections, the Word report writers (`officer`, `flextable`), the atlas tools (`magick`, `pdftools`) and the dependency installer itself (`BiocManager`). Missing Bioconductor packages are installed through `BiocManager` and CRAN packages through `install.packages()`, each failing package is retried alone so its build log can be reported, and the versions actually used are written to `provenance/packages.csv`. `Rscript scripts/install_packages.R` runs the same installation without starting an analysis.
- The run records the configuration, input checksums, sample metadata, contrast summary, and session information under `provenance/`.

## Output layout

Each run writes to `results/<name of the data file>/`, so two experiments never share a folder:

```text
results/NIH3T3/          # data file: NIH3T3.xlsx
├── 00_differential_expression/   DESeq2 tables, normalized counts, QC figures
├── section01_screening/          DEG tables, sData.xlsx, volcano/FC-Exp panels, heat maps, chromosome landscape
├── section02_go_ora/             GO ORA tables, dot plots, rich-factor plot, functional-class map, ten per-class term-gene networks
├── section03_gsea_go/            GSEA tables, per-set PDFs, composite panels, enrichment maps, MSigDB Hallmark
├── section04_kegg_pathview/      KEGG table, source note, lollipop/UpSet panels, pathway maps and KGML
├── section05_tf_regulon/         TF edges, regulon summary (with tf_fpkm / expressed), TF dot plot, TF expression dot plot, chromosome maps of the top-10 regulons
├── section06_regulatory_activity/ TF activity (ULM/AUCell/WMean), GSVA pathway activity
├── section07_tf_network/         TF–target network figures and node/edge tables
├── section08_tf_modules/         TF module network, membership/similarity/module tables
├── section09_geneset_heatmaps/   GO heat maps and TF-regulon dot plots (≤30 genes each)
├── section10_ppi_network/        STRING PPI network, hub-protein ranking, module GO annotation, edge/node/module tables
├── section11_subcellular_localization/  subcellular localization map, compartment and gene tables
├── provenance/                   run manifest, sessionInfo, and config_used.yml (the resolved settings)
└── report/                       rendered English manuscript draft (self-contained HTML),
                                  fixed-template English Word report (results_report.docx) and its values table
```

The DoRothEA collections live in `data/dorothea/` (`dorothea_hs.rda`, `dorothea_mm.rda`), mirroring the objects exported by the `dorothea` package. The demonstration inputs are `input/NIH3T3.xlsx` and `input/NIH3T3_metadata.csv`; no other file is needed for a run.

Rebuild the figure atlas with:

```bash
Rscript scripts/build_figure_atlas.R                 # newest run under results/ (writes RNAseq_results_atlas.pdf)
Rscript scripts/build_figure_atlas.R results/my_run  # a specific one
```

## Verifying the run

The shipped demonstration was executed end to end on this host with R 4.6.1 / Bioconductor 3.23 (`DESeq2` 1.52.0, `clusterProfiler` 4.20.0, `msigdbr` 26.1.1, `pathview` 1.52.0, `apeglm` 1.34.0, `decoupleR` 2.17.0, `GSVA` 2.6.6, `org.Mm.eg.db` 3.23.0), through `run_all(counts_file = "input/NIH3T3.xlsx", species = "mouse", log2fc = 0.75, min_fpkm = 5, top_tfs = 40)`. Every number and figure under `results/NIH3T3/` is computed output for the six supplied samples, not a placeholder. The sample call above uses `top_tfs = 30`, which changes only how many transcription factors are drawn; the differential-expression and enrichment results are unchanged.

`Rscript scripts/install_packages.R` installs any missing package from the same list the run checks (Bioconductor through `BiocManager`, CRAN through `install.packages`); it reports the build log of anything that still fails.

## Raw reads and publication readiness

The supplied workbook begins at gene-level counts. For a complete raw-read workflow, archive FastQC/MultiQC, trimming, alignment or pseudoalignment, reference genome and annotation versions, count construction, mapping metrics, and the exact workflow version. The upstream specification is in [docs/raw_reads_to_counts.md](docs/raw_reads_to_counts.md).

Before submission, add raw-read accession numbers, source data, a stable code repository, biological validation, and a fully reviewed statistical design. The English manuscript template is in [manuscript/manuscript_draft.Rmd](manuscript/manuscript_draft.Rmd).

## Reproducing the paper

The repository contains everything needed to rebuild the case-study paper from
the released inputs:

```bash
# 1. the analysis (eleven sections, figure atlas, Word report)
Rscript pipeline/run_analysis.R --counts input/NIH3T3.xlsx \
        --samples input/NIH3T3_metadata.csv --species mouse \
        --contrast FBS,Control

# 2. the tables, taken directly from results/NIH3T3/
python3 work/make_paper_tables.py --results results/NIH3T3

# 3. the composite figures, the Word build and the PDF build
python3 work/build_word_figures.py
python3 work/make_docx_md.py && python3 work/postprocess_word.py
python3 work/build_manuscript_pdf.py
```

`work/` holds those build scripts and their own [README](work/README.md); the
analysis itself never calls them. Every table in the paper is regenerated from
the result tables of the run between the `<!--TABLE:n-->` markers, and Table 8
is a literature comparison that is maintained by hand.

## Citation

If you use TransAtlas, please cite the software and the paper:

```bibtex
@software{TransAtlas,
  author    = {Sun, Changye},
  title     = {TransAtlas: an Integrated R Workflow for Downstream Bulk RNA-seq Analysis},
  year      = {2026},
  version   = {1.2.0},
  url       = {https://github.com/hscsun/TransAtlas-R-v1.2}
}
```


## Licence

Released under the [MIT licence](LICENSE) © 2026 Changye Sun. The workflow calls
third-party R packages and data resources that carry their own licences
(Bioconductor packages, DoRothEA, MSigDB, KEGG, STRING, the Gene Ontology);
users are responsible for observing them when redistributing results.

## Contact

Changye Sun — Henan Key Laboratory of Tissue Regeneration, Henan Medical
University, Xinxiang, Henan 453003, China · sunchangye-boy@163.com.
Bugs and feature requests: [open an issue](https://github.com/ChangyeSun/TransAtlas-R/issues).
