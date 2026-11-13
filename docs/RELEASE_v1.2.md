# TransAtlas-R v1.2

Prepared by ClawsGO Science Agent. Release date: 2026-11-12.

The distribution is named `TransAtlas-R-v1.2.zip`; its semantic software version is **1.2.0**. VERSION, README, CITATION.cff, Zenodo metadata and the run manifest use 1.2.0. This is a local release distribution; no GitHub or Zenodo publication was performed for this update.

## Unbounded calculation effects and independent cutoffs

All calculations, gene rankings and effect columns retain true unshrunk log2FC without an upper bound. This includes screening CSV/XLSX tables, GO/MSigDB/KEGG/TF GSEA inputs, TF-target edges, mean/median target effects, module summaries and network-node exports. Shrinkage estimates remain separate inspection tables.

| Parameter | Default | Scope |
|---|---:|---|
| `up_cutoff` | 2.5 | Gene displays only: symmetric coordinates, chromosome stems and colour scales from -2.5 to 2.5. NULL uses an automatic range. |
| `low_cutoff` | 0.15 | GSEA input: absolute true log2FC >= 0.15, together with gene BH-adjusted P < alpha. NULL or 0 removes the effect-size lower bound. |
| `log2fc` | 0.75 | Primary DEG calls and DEG gene sets. Independent of both cutoff parameters above. |

Changing the display limit or screening ladder leaves the GSEA ranking unchanged. `log2fc_cap` is accepted as a deprecated display-only alias. Command-line forms are `--up-cutoff` and `--low-cutoff`; use `--low-cutoff 0` to remove the lower bound. The low cutoff is a minimum absolute effect, so the default excludes effects strictly between -0.15 and 0.15 while retaining both tails.

`provenance/gsea_input_genes.csv`, `gsea_ranked_vector.csv` and `status_gsea_input.json` record the actual input. ORA backgrounds use all eligible tested genes with usable adjusted P and identifiers after expression filtering, independently of GSEA significance/effect screening.

Related fixes found in the same review:

- Duplicate identifiers retain the finite effect with largest absolute magnitude, including negative effects.
- TF module counts and mean effects count each target once per TF.
- The TF odds ratio uses the full 2x2 contingency table.
- Chromosome-map thresholds are included exactly even with a custom screening ladder; automatic display ranges follow true effects.
- The optional legacy expression filter is applied consistently to the eligible universe and DEG sets.
- Missing network effects remain NA; node-size captions match the fixed display sizes.
- The standalone GO network defaults to unshrunk effects, accepts `--up-cutoff`, and uses graphlayouts for its stress layout.
- Reports keep nonsignificant all-term GO fallback rows out of significance claims, use the correct GO table path and describe the configured primary heat-map threshold. The localization narrative uses the annotated DEG denominator instead of the largest compartment count; key-target references and separate shrinkage-table descriptions are corrected.

These changes affect GSEA enrichment and ORA/regulatory results. Counts and estimates from the earlier capped-ranking distribution are historical comparisons, not the current reference.

## Other v1.2 changes

### Activity profiles remain auditable

Differential TF and GSVA activity tables retain every scored profile. The existing mean, difference, statistic, P-value and adjusted-P columns are preserved. New columns are `n_reference`, `n_treatment`, `test_status` and `test_reason`.

- `tested`: a finite two-sided Welch t-test P value is available.
- `not_testable`: fewer than two finite observations in either group, zero/numerically negligible within-group variation, or a non-finite test result. Group means and differences remain available; P values and adjusted P values are NA.
- `error`: an unexpected test error, with its message retained.

BH adjustment includes only successfully tested profiles. No synthetic P values or replacement significance tests are introduced. Figures requiring P values use valid tests; the full differential table retains the other profiles. Report text identifies the profiles without a valid test.

For NIH3T3, all three TF methods retain 87 profiles. ULM and WMean test all 87. AUCell tests 76; eleven have constant within-group scores. Nfic previously disappeared after a t-test error; ten additional constant profiles previously had NaN P values. All eleven now have explicit status and reason. P values for the 76 previously valid AUCell tests are unchanged, and 52 AUCell TFs remain significant at adjusted P < 0.05.

### Exact MSigDB source

Collection availability is checked before retrieval. `H` retains the human Hallmark collection and msigdbr ortholog mapping for mouse, maintaining the previous analysis policy. `MH` explicitly selects mouse-native Hallmark when supported by the installed database. H and MH remain distinct choices.

Section 3 writes `GSEA_MSigDB_<collection>_source.json`; section 6 writes `GSVA_<collection>_source.json`. They record requested and actual databases, database and package versions, collection/subcollection, database and target species, ortholog mapping, set/gene counts, fallback attempts and an MD5 checksum of the sorted unique set–symbol membership. The checksum identifies the membership used; it is not a cryptographic security assertion.

For the validated default mouse run, both use MSigDB **2026.1.Hs**, collection **H**, mapped to **Mus musculus** through msigdbr 26.1.1: 50 sets, 4,393 unique symbols and 7,348 membership pairs. This corrects the earlier broad statement that Hallmark lacks a mouse-native collection: the installed mouse database does provide **MH**.

### Run warning records

Each invocation creates a fresh `warnings.csv` in the output root with run ID, sequence, UTC time, stage, condition class, message and condition call where supplied by R. Each visited stage folder also receives a warning table; stages sharing a folder preserve all records from that invocation. Reports feed their intercepted warnings into the same central log. A warning summary and `status_warnings.json` are finalized on successful completion or an R error. An abrupt external process termination may prevent finalization; warnings already appended remain on disk.

Warnings retain their normal R behavior. The recorder observes warnings reaching the pipeline handler and explicitly integrates report warnings; operations deliberately suppressing or internally consuming warnings can still limit what reaches it. Optional errors remain in `optional_errors.log`. The warning log is independent of failure status: a completed run can have documented warnings.

### Drawing fixes found during warning review

- Apply legend size overrides once when colour/fill/shape guides merge.
- Explicitly allow the intentionally unlabelled GO hub in the text layer.
- Return NA chromosome count when no mapping is available, avoiding an empty maximum.

These changes preserve the plotted data and statistical inputs.

## Validation

The completed rerun is `results/NIH3T3_v1.2_uncapped`, with `up_cutoff = 2.5`, `low_cutoff = 0.15`, `log2fc = 0.75`, `min_fpkm = 5` and `top_tfs = 30`. All eleven sections, QC, the overall run and both automatic reports completed. Validation passed 25 targeted cutoff checks, two CLI checks, five report edge fixtures and 43 full-result checks. An actual standalone stress-layout network completed with a display limit of 1.2 while its exported effects exceeded that limit. The Word report was rendered to 21 PDF pages and visually reviewed. The final report was rebuilt from the completed tables after correcting its localization denominator and wording; this step is recorded in `uncapped_report_rebuild.log`. Checks, logs and the warning assessment are included in the separate validation bundle.

The GSEA vector contains 5,412 unique identifiers with true log2FC ranging from **-8.45911971309203 to 2.76626656240079**. Default gene padj and lower-cutoff requirements are checked directly. Changing the display cutoff or screening ladder does not change the vector; low_cutoff NULL/0 admits additional small significant effects while retaining the padj filter. The primary DEG call remains **1,213 genes (830 up, 383 down)**. Current GO ORA counts are **184 / 42 / 63** for BP/MF/CC; GO GSEA counts are **405 / 61 / 113**; Hallmark has **16 of 46 tested** significant sets and KEGG has **27** significant pathways.

TF analysis tested 63 expressed regulons: 26 were significant by ORA and 16 by GSEA. The complete table contains 2,304 TF-target edges; all target effects and all 63 regulon mean/median effects were checked against eligible unshrunk estimates. The STRING network retains 236 proteins, 1,086 interactions and ten modules. Localization evaluated 1,100 CC-annotated DEGs across 15 compartments; 1,057 of these genes map into at least one of the selected high-level compartments.

The run recorded **49 warnings**, all assessed below. No optional errors, capped-ranking tie warnings or previously fixed drawing warnings were recorded.

| Warning category | Occurrences | Interpretation |
|---|---:|---|
| Annotation coverage | 39 | Unmatched symbols have no OrgDb mapping; no substitute identifiers or effects are invented. Repeated per-regulon background mapping accounts for most occurrences. |
| Automatic semantic-data setup | 3 | The GO dependency computed semantic data automatically; simplification completed. |
| Very small enrichment P estimates | 5 | Estimates reached the configured backend numerical floor. Interpret these values at that limit. |
| ggupset dependency deprecation | 2 | The dependency uses an older line-size aesthetic; both plots completed. |

Pathview remote map downloads are disabled for validation. Available packages, a valid DESeq2 fit and KEGG/STRING resource caches are reused. Clean dependency installation, fresh remote downloads, other operating systems and human/rat end-to-end runs are outside this rerun. Numerical tables and composites of the static submission paper remain the earlier case-study snapshot and carry a dated note. The automated HTML and Word reports use the current results.

## Running

The README Sample run and public entry point are unchanged:

```r
source("pipeline/run_all.R")
run_all(
  counts_file = "input/NIH3T3.xlsx",
  species = "mouse",
  log2fc = 0.75,
  min_fpkm = 5,
  up_cutoff = 2.5,
  low_cutoff = 0.15,
  top_tfs = 30
)
```

The source distribution includes previous review fixes, the NIH3T3 example and bundled DoRothEA resources. The validation bundle is separate from production sources and excludes large fit caches.
