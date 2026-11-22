"""Package existing scalar results. Requires pandas, NumPy, openpyxl and PyMuPDF."""
import argparse
import csv
from collections import Counter
from decimal import Decimal
from hashlib import sha256
from io import BytesIO
import math
from pathlib import Path
import re
from zipfile import ZipFile

import fitz
import numpy as np
import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT
NAME = "supplementary_numerical_results.xlsx"
OBS = "observed_data/"
REC = "semi_synthetic/"
SENS = REC
FIG = "semi_synthetic/"
RUN = "diagnostics/"
OM, OF, RM, RF, PAIR, RT = (
    "OBSERVED_METRICS", "OBSERVED_FIT_COUNTS", "RECOVERY_METRICS",
    "RECOVERY_FIT_COUNTS", "PAIRED_METRIC_COUNTS", "RUNTIME",
)
OBS_KEYS = ["dataset", "model", "style"]
REC_KEYS = ["dataset", "model", "dropped_modality", "mechanism", "target_rate"]
OBS_LABELS = {
    "rand_diag": "Diagnosis / disease-group Rand index",
    "ari_diag": "Diagnosis / disease-group adjusted Rand index",
    "rand_miss": "Modality-availability Rand index",
    "ari_miss": "Modality-availability adjusted Rand index",
    "chi_diag": "Diagnosis / disease-group Calinski–Harabasz index",
    "db_diag": "Diagnosis / disease-group Davies–Bouldin index",
    "chi_miss": "Modality-availability Calinski–Harabasz index",
    "db_miss": "Modality-availability Davies–Bouldin index",
}
REC_LABELS = {
    "rv": "Latent-score RV", "rv_W": "Loading RV",
    "gw": "Gromov–Wasserstein distance", "wgene_rankcor": "Loading-rank Spearman correlation",
    "wgene_jaccard": "Top-50 Jaccard overlap", "recon_cor": "Reconstruction correlation",
    "recon_relerr_full": "Relative Frobenius error (overall)",
}
MODALITIES = {"gene": "Gene expression", "meta": "Metabolomics", "mrna": "mRNA",
              "methylation": "DNA methylation", "rppa": "RPPA"}
REC_LABELS.update({f"recon_relerr_{k}": f"Relative Frobenius error ({v})"
                   for k, v in MODALITIES.items()})
METHODS = dict(zip(
    "mofa cmf rgcca sgcca gpcca block_pls gipca bidifac bsf nemo icluster fabia gfa "
    "o2pls cia mcia jive slide tica intnmf".split(),
    "MOFA CMF RGCCA sGCCA GPCCA Block-PLS GIPCA BIDIFAC BSF NEMO iCluster FABIA GFA "
    "O2PLS CIA MCIA JIVE SLIDE tICA intNMF".split(),
))
ID = ["Dataset", "Method"]
CONFIG = ID + ["Removed modality", "Missingness mechanism", "Target missingness rate"]
COUNTS = ["Number attempted", "Number successful", "Number failed", "Success rate"]


def number(value):
    return Decimal(value) if value and value.lower() not in {"na", "nan", "inf", "-inf"} else None


def integer(value):
    parsed = number(value)
    return int(parsed) if parsed is not None else None


def rate(value):
    # Raw CSVs serialize the same two target levels with different binary round-trip spellings.
    parsed = Decimal(value)
    target = min((Decimal(".2"), Decimal(".7")), key=lambda item: abs(item - parsed))
    if abs(target - parsed) > Decimal("1e-12"):
        raise ValueError(f"Unexpected target rate {value}")
    return target


def key(row, columns):
    return tuple(rate(row[c]) if c == "target_rate" else
                 row[c].lower() if c == "dataset" else row[c] for c in columns)


def identity(row):
    return [row["dataset"].upper(), METHODS[row["model"]]]


def configuration(row):
    return identity(row) + [MODALITIES[row["dropped_modality"]],
                            row["mechanism"], rate(row["target_rate"])]


def finite(series):
    return series[np.isfinite(series)]


def fraction(success, attempted):
    return Decimal(repr(success / attempted)) if attempted else None


def main(destination, manuscript=None, supporting_pdf=None):
    targets = [destination / NAME, destination / "source_manifest.csv",
               destination / "verification.txt"]
    if any(path.exists() for path in targets):
        raise FileExistsError("Choose an output directory without existing deliverables.")
    provenance_path = PACKAGE / "docs/result_provenance.csv"
    with provenance_path.open(newline="") as handle:
        provenance = {r["file"]: r for r in csv.DictReader(handle)}
    manifest = []
    tables = {}
    checks = Counter()
    discrepancies = []

    def source(sheet, relative, role):
        relative = "results/" + relative + ".csv"
        path = PACKAGE / relative
        digest = sha256(path.read_bytes()).hexdigest()
        meta = provenance.get(relative, {})
        if meta.get("file_sha256") and digest != meta["file_sha256"]:
            discrepancies.append(f"Provenance checksum differs for {relative}.")
        manifest.append([sheet, role, str(path.relative_to(ROOT)), digest,
                         meta.get("source_file", ""), meta.get("source_sha256", "")])
        return path

    def records(sheet, relative, role="Export values"):
        with source(sheet, relative, role).open(newline="") as handle:
            return list(csv.DictReader(handle))

    def frame(sheet, relative, role="Recalculate from saved scalar records"):
        data = pd.read_csv(source(sheet, relative, role), keep_default_na=False,
                           na_values=["", "NaN", "Inf", "-Inf"], float_precision="round_trip",
                           dtype={"target_rate": str})
        if "target_rate" in data:
            data["target_rate"] = data["target_rate"].map(rate)
        return data

    def compare(label, actual, expected, location=""):
        checks[label] += 1
        if actual is None or expected is None:
            equal = actual is expected
        else:
            equal = math.isclose(float(actual), float(expected), rel_tol=2e-12, abs_tol=2e-13)
        if not equal:
            discrepancies.append(f"{label} {location}: exported {actual}; reference {expected}.")

    observed = records(OM, OBS + "bootstrap_metric_summary")
    obs_raw = frame(OM, OBS + "bootstrap_metrics_by_replicate")
    assert not obs_raw.duplicated(OBS_KEYS + ["boot"]).any()
    obs_groups = {idx: group for idx, group in obs_raw.groupby(OBS_KEYS)}
    obs_index = {}
    rows = []
    for row in observed:
        metric = row["metric"]
        if metric not in OBS_LABELS and metric != "R2":
            continue
        idx = key(row, OBS_KEYS)
        good = obs_groups[idx].query("status in ['success', 'partial'] and n_successful > 0")
        values = finite(good[metric])
        if row["style"] == "CC_SF" and metric.endswith("_miss"):
            values = values.iloc[:0]
        n = integer(row["n"])
        estimates = [number(row[c]) for c in ("mean", "low", "high")]
        compare("Observed finite count", n, len(values), (idx, metric))
        expected = [values.mean(), *values.quantile([.025, .975])] if len(values) else [None] * 3
        for actual, reference in zip(estimates, expected):
            compare("Observed bootstrap statistic", actual, reference, (idx, metric))
        obs_index[(*idx, metric)] = (n, *estimates)
        label = (("MMSE" if idx[0] == "adni" else "Karnofsky") + " R²"
                 if metric == "R2" else OBS_LABELS[metric])
        status = ("Undefined for complete cases" if row["style"] == "CC_SF" and metric.endswith("_miss")
                  else "No finite eligible replicates" if not n else None)
        rows.append(identity(row) + [row["style"], label, *estimates, n, status])
    tables[OM] = (ID + ["Missing-data strategy", "Evaluation metric", "Mean metric value",
                        "Bootstrap 2.5th percentile", "Bootstrap 97.5th percentile",
                        "Valid bootstrap replicates", "Availability"], rows)
    for file, columns, label in [
        ("observed_data/clinical_alignment_bootstrap_summary", ("mean", "low", "high"), "Figures 2–3"),
        (OBS + "adjusted_rand_plotted_values", ("mean", "p025", "p975"), "S5"),
    ]:
        for row in records(OM, file, "Figure cross-check"):
            idx = key(row, OBS_KEYS) + (row["metric"],)
            if idx in obs_index:
                for actual, expected in zip(obs_index[idx], [integer(row["n"]), *[number(row[c]) for c in columns]]):
                    compare(label, actual, expected, idx)

    reliability = records(OF, "diagnostics/reliability_all_methods_datasets_strategies")
    saved_counts = {key(r, OBS_KEYS): r for r in records(
        OF, OBS + "bootstrap_fit_success_counts", "Count cross-check")}
    figure_counts = {key(r, OBS_KEYS): r for r in records(
        OF, "diagnostics/fitting_reliability_plotted_values", "Figure S6 cross-check")}
    rows = []
    for row in reliability:
        idx = key(row, OBS_KEYS)
        supported = row["supported"] == "True"
        group = obs_groups.get(idx)
        attempted, successful, failed = [integer(row[f"recorded_{c}_units"])
                                         for c in ("attempted", "successful", "failed")]
        absent = integer(row["missing_historical_records"])
        prefit = integer(row["pre_fit_failed_slots"])
        if supported:
            saved = saved_counts[idx]
            figure = figure_counts[idx]
            for c, actual in zip(("attempted", "successful", "failed"), (attempted, successful, failed)):
                compare("Observed counters vs summary", actual, integer(saved[c + "_fits"]), idx)
                compare("Observed counters vs saved records", actual, group["n_" + c].sum(), idx)
            compare("S2 attempts", attempted, integer(figure["attempted_fits"]), idx)
            compare("S2 successes", successful, integer(figure["successful_fits"]), idx)
            compare("Observed historical absence", absent, 1000 - len(group), idx)
        level = ("Bootstrap summary" if row["counter_unit"] == "normalized bootstrap summary"
                 else "Individual fit / workflow slot" if supported else None)
        status = ("Strategy absent from saved workflow" if not supported
                  else "One pre-fit screening failure recorded as five slots" if prefit
                  else "Recorded failures or partial records" if failed else "All recorded units successful")
        rows.append(identity(row) + [row["style"], attempted if supported else None,
            successful if supported else None, failed if supported else None,
            number(row["recorded_success_rate"]), level,
            *[integer(row[c]) if supported else None for c in (
                "recorded_bootstrap_records", "fully_successful_records", "partial_records", "failed_records",
                "missing_historical_records", "pre_fit_screening_failed_records", "pre_fit_failed_slots")], status])
    tables[OF] = (ID + ["Strategy"] + COUNTS + ["Count level", "Recorded bootstrap summaries",
                  "Fully successful summaries", "Partial summaries", "Failed summaries",
                  "Missing historical summaries", "Pre-fit screening failures", "Pre-fit failed slots",
                  "Failure / status information"], rows)
    compare("Manuscript Methods missing summaries", sum(r[12] or 0 for r in rows), 48)

    primary = frame(RM, REC + "recovery_metrics_by_repetition")
    sensitivity = frame(RM, SENS + "intnmf_sensitivity_metrics_by_repetition")
    assert not primary.duplicated(REC_KEYS + ["repetition"]).any()
    assert not sensitivity.duplicated(REC_KEYS + ["repetition"]).any()
    main_groups = {idx: g for idx, g in primary.groupby(REC_KEYS)}
    sensitivity_groups = {idx: g for idx, g in sensitivity.groupby(REC_KEYS)}
    rows = []
    recovery_index = {}
    for file, groups in [(REC + "recovery_metric_summary", main_groups),
                         (SENS + "intnmf_sensitivity_metric_summary", sensitivity_groups)]:
        for row in records(RM, file):
            idx = key(row, REC_KEYS)
            if row["metric"].startswith("recon_relerr_") and row["metric"] != "recon_relerr_full":
                applicable = {"gene", "meta"} if idx[0] == "adni" else {"mrna", "methylation", "rppa"}
                if row["metric"].removeprefix("recon_relerr_") not in applicable:
                    assert not integer(row["n"])
                    continue
            group = groups[idx]
            metric = row["metric"]
            values = finite(group.loc[group.status.eq("success"), metric])
            n = integer(row["n"])
            mean = number(row["mean"])
            reference = values.mean() if len(values) else None
            compare("Recovery metric finite count", n, len(values), (idx, metric))
            compare("Recovery metric mean", mean, reference, (idx, metric))
            sd = Decimal(repr(float(values.std(ddof=1)))) if len(values) > 1 else None
            experiment = row.get("experiment", "Main benchmark")
            recovery_index[(*idx, metric)] = (n, mean)
            rows.append(configuration(row) + [experiment, REC_LABELS[metric], mean, sd, n])
    tables[RM] = (CONFIG + ["Experiment", "Evaluation metric", "Mean", "Sample standard deviation",
                           "Valid repetitions"], rows)
    for row in records(RM, "semi_synthetic/intnmf_recovery_plotted_values", "S10 cross-check"):
        idx = key(row, REC_KEYS) + (row["metric"],)
        n, mean = recovery_index[idx]
        compare("S10 finite count", n, integer(row["successful_metric_repetitions"]), idx)
        compare("S10 mean", mean, number(row["value"]), idx)
    for row in records(RM, FIG + "feature_recovery_plotted_values", "S6 cross-check"):
        means = []
        for mechanism in ("MCAR", "MAR", "MNAR"):
            idx = (row["dataset"], row["model"], row["modality"], mechanism, rate(row["rate"]), row["metric"])
            n, mean = recovery_index[idx]
            compare("S6 mechanism finite count", n, integer(row["n_" + mechanism]), idx)
            compare("S6 mechanism mean", mean, number(row["mean_" + mechanism]), idx)
            if mean is not None:
                means.append(mean)
        compare("S6 equally weighted mean", sum(means) / len(means) if means else None, number(row["value"]), idx)

    source(RF, REC + "recovery_metrics_by_repetition", "Failure stages and counter cross-check")
    source(RF, SENS + "intnmf_sensitivity_metrics_by_repetition", "Sensitivity counters and failure stages")
    recovery_counts = {key(r, REC_KEYS): r for r in records(RF, REC + "recovery_fit_success_counts")}
    baseline_frames = [
        frame(RF, REC + "tcga_cia_o2pls_baseline_failures", "Excluded baseline status records"),
        frame(RF, SENS + "tcga_intnmf_original_baseline_failure_records", "Original excluded intNMF baseline status records"),
    ]
    baseline = pd.concat(baseline_frames, ignore_index=True)
    assert set(baseline.status) == {"baseline_failed"}
    rows = []
    for groups in (main_groups, sensitivity_groups, dict(tuple(baseline.groupby(REC_KEYS)))):
        for idx, group in groups.items():
            statuses = group.status.value_counts()
            success = int(statuses.get("success", 0))
            stages = [int(statuses.get(s, 0)) for s in
                      ("imputation_failed", "fit_failed", "metric_failed", "baseline_failed")]
            assert len(group) == success + sum(stages)
            slots = int(group.n_attempted.sum())
            attempted = slots - stages[0] - stages[3]
            failed = stages[1] + stages[2]
            assert attempted == success + failed
            if groups is main_groups:
                saved = recovery_counts[idx]
                for c, value in [("records", len(group)), ("successful", success), ("failed", failed)]:
                    compare("Recovery count summary", value, integer(saved[c]), idx)
            row = dict(zip(REC_KEYS, idx))
            experiment = ("Original TCGA baseline (excluded)" if stages[3]
                          else group.experiment.iloc[0] if groups is sensitivity_groups else "Main benchmark")
            rows.append(configuration(row) + [experiment, attempted, success, failed, fraction(success, attempted),
                         len(group), *stages])
    tables[RF] = (CONFIG + ["Experiment"] + COUNTS + ["Recorded repetition slots",
                  "Imputation failures", "Refit failures", "Metric-calculation failures",
                  "Baseline failure slots"], rows)
    for dataset, methods, success, patients in [("adni", 13, 15562, 450), ("tcga", 11, 19791, 465)]:
        subset = primary[primary.dataset.eq(dataset)]
        compare("Manuscript p.23 primary methods", subset.model.nunique(), methods, dataset)
        compare("Primary recovery successes", subset.status.eq("success").sum(), success, dataset)
        compare("Manuscript Methods recovery sample size", subset.n_samples.min(), patients, dataset)
        compare("Manuscript Methods recovery sample size", subset.n_samples.max(), patients, dataset)

    source(PAIR, REC + "recovery_metrics_by_repetition", "Verify paired finite means and counts")
    paired_keys = set()
    rows = []
    for row in records(PAIR, FIG + "patient_feature_recovery_plotted_values"):
        row.update(dropped_modality=row["modality"], target_rate=row["rate"])
        idx = key(row, REC_KEYS)
        paired_keys.add(idx)
        group = main_groups[idx]
        paired = group.loc[group.status.eq("success") & np.isfinite(group.rv) & np.isfinite(group.wgene_jaccard)]
        n, rv, jaccard = integer(row["n"]), number(row["rv"]), number(row["jaccard"])
        compare("S7 paired count", n, len(paired), idx)
        compare("S7 paired latent RV mean", rv, paired.rv.mean(), idx)
        compare("S7 paired Jaccard mean", jaccard, paired.wgene_jaccard.mean(), idx)
        rows.append(configuration(row) + [rv, jaccard, n])
    eligible = {idx for idx, group in main_groups.items()
                if (group.status.eq("success") & np.isfinite(group.rv) & np.isfinite(group.wgene_jaccard)).any()}
    assert paired_keys == eligible and len(rows) == len(paired_keys)
    for ds, expected in [("ADNI", 78), ("TCGA", 99)]:
        for target in (Decimal(".2"), Decimal(".7")):
            compare("S7 panel configuration count", sum(r[0] == ds and r[4] == target for r in rows), expected, (ds, target))
    tables[PAIR] = (CONFIG + ["Mean latent-score RV", "Mean top-50 Jaccard overlap",
                             "Repetitions with both metrics finite"], rows)

    def runtime_stats(row, values, idx):
        values = finite(values)
        values = values[values.gt(0)]
        n = integer(row["n_timed"])
        quantiles = [number(row[c]) for c in ("median_seconds", "q25_seconds", "q75_seconds")]
        compare("Runtime timed count", n, len(values) if n is not None else None, idx)
        expected = list(values.quantile([.5, .25, .75])) if len(values) else [None] * 3
        for a, b in zip(quantiles, expected):
            compare("Runtime quantile", a, b, idx)
        return [*quantiles, n]

    source(RT, OBS + "bootstrap_metrics_by_replicate", "Verify elapsed-time eligibility and quantiles")
    source(RT, REC + "recovery_metrics_by_repetition", "Verify elapsed-time eligibility and quantiles")
    rows = []
    runtime_index = {}
    for row in records(RT, RUN + "observed_runtime_summary"):
        row["style"] = row["strategy"]
        idx = key(row, OBS_KEYS)
        group = obs_groups.get(idx)
        fits = integer(row["fits_per_record"])
        values = pd.Series(dtype=float)
        if group is not None:
            eligible = group.status.eq("success") & group.n_attempted.eq(fits) & group.n_successful.eq(fits) & group.n_failed.eq(0)
            values = group.loc[eligible, "elapsed_seconds"] / fits
        stats = runtime_stats(row, values, idx)
        runtime_index[idx] = stats
        rows.append(["Observed", *identity(row), row["style"], None, None, None,
                     *stats, fits, integer(row["n_recorded"]), row["coverage"]])
    pooled_keys = [c for c in REC_KEYS if c != "mechanism"]
    pooled = {idx: g for idx, g in primary.groupby(pooled_keys)}
    for row in records(RT, RUN + "semisynthetic_runtime_summary"):
        idx = key(row, pooled_keys)
        group = pooled[idx]
        stats = runtime_stats(row, group.loc[group.status.eq("success"), "fit_elapsed_seconds"], idx)
        rows.append(["Semi-synthetic", *identity(row), None, MODALITIES[row["dropped_modality"]],
                     "MCAR / MAR / MNAR pooled", rate(row["target_rate"]), *stats, 1,
                     integer(row["n_recorded"]), "Available"])
    tables[RT] = (["Analysis type"] + ID + ["Strategy", "Removed modality", "Missingness mechanism",
                  "Target missingness rate", "Median runtime (seconds)", "First quartile (seconds)",
                  "Third quartile (seconds)", "Successful timed records used",
                  "Fits per timed record", "Recorded summaries or repetition slots", "Availability"], rows)

    if manuscript is not None:
        with fitz.open(manuscript) as pdf:
            pages = [" ".join(page.get_text().split()) for page in pdf]
        assert len(pages) == 37
        assert "Figure 8: Observed-data runtime" in pages[13]
        assert "472 of 1,000" in pages[7] and "639 of 1,000" in pages[7]
    if supporting_pdf is not None:
        with fitz.open(supporting_pdf) as pdf:
            assert len(pdf) in (10, 37)
            page = 8 if len(pdf) == 10 else 35
            assert "S9 Fig. Semi-synthetic refitting runtime" in " ".join(pdf[page].get_text().split())
    compare("Manuscript observed results ADNI GIPCA NA successes", saved_counts["adni", "gipca", "NA"]["successful_fits"], 472)
    compare("Manuscript observed results TCGA BIDIFAC NA successes", saved_counts["tcga", "bidifac", "NA"]["successful_fits"], 639)
    compare("Manuscript observed results TCGA GPCCA CC failures", saved_counts["tcga", "gpcca", "CC_SF"]["failed_fits"], 1000)
    compare("Manuscript observed results ADNI GPCCA runtime (rounded)", round(runtime_index["adni", "gpcca", "NA"][0]), 6638)
    compare("Manuscript observed results ADNI MOFA runtime (rounded)", round(runtime_index["adni", "mofa", "NA"][0], 1), 7.2)

    readme = [
        ("Scope", "Saved final benchmarking results for ADNI and TCGA. This workbook summarizes the observed bootstrap and semi-synthetic benchmarks. See source_manifest.csv for worksheet-level CSV provenance."),
        ("Missing values", "Blank estimate cells indicate an undefined or unavailable quantity. A zero valid count indicates that no eligible finite values were saved. Unsupported observed strategies have blank fit counts."),
        ("Numeric precision", "Source decimals are stored at full precision in the workbook XML. Excel displays and calculates with its own numeric precision (about 15 significant digits). Missingness rates are the nominal design levels, displayed as percentages."),
        ("Observed strategies", "CC_SF uses complete cases and screened features. MI_SF uses multiple imputation and screened features. NA fits directly to incomplete data."),
        ("Observed metrics", "Rand index is pairwise agreement. Adjusted Rand corrects for chance. Both compare latent-space clusters with diagnosis/disease group or modality-availability patterns. Availability metrics are undefined for CC_SF."),
        ("Observed metrics", "Calinski–Harabasz is the between- to within-group dispersion ratio with degrees-of-freedom adjustment. Davies–Bouldin measures average worst-case cluster similarity based on dispersion and separation. MMSE (ADNI) and Karnofsky (TCGA) R² describe the in-sample clinical outcome regression on latent scores."),
        ("Observed denominators", "Means and empirical 2.5th/97.5th percentiles use finite values from recorded successful or partially successful bootstrap summaries with at least one successful fit. MI_SF metrics are saved within-bootstrap averages across successful imputations. Percentiles use linear interpolation."),
        ("Observed reliability", "R and ADNI Python MI_SF use individual-fit/workflow-slot counters. Each ADNI MOFA and CMF screening failure accounts for five failed workflow slots before model fitting. TCGA MOFA and CMF counters are normalized bootstrap summaries for all strategies. Individual imputation-fit counts are unavailable there."),
        ("Observed historical gaps", "Of the planned observed summaries, 48 are missing in ADNI: five per strategy for MOFA and CMF and one MI_SF summary for each of 18 R methods. Missing summaries are excluded from attempt and failure denominators."),
        ("Observed failure", "TCGA GPCCA CC_SF recorded 1,000 failures and zero successes. The manuscript links initialization failure to singular covariance involving duplicate methylation features."),
        ("Recovery design", "MCAR is missing completely at random, MAR is missing at random conditional on retained assays, and MNAR is missing not at random conditional on the removed assay. Rates are target probabilities. Each main configuration has 100 recorded repetition slots."),
        ("Recovery metrics", "Latent-score RV and loading RV measure matrix association on comparable scores or loadings. Gromov–Wasserstein distance compares latent geometries. Spearman compares feature importance ranks; Jaccard measures intersection over union for the top 50 features. Importance is based on loading norms."),
        ("Reconstruction", "Reconstruction correlation is the overall vectorized Pearson correlation. Relative Frobenius error is the error norm divided by the reference norm, reported overall and by modality. The saved workflow reports reconstruction correlation at the overall level."),
        ("Recovery denominators", "Means and sample standard deviations use the successful finite values for each metric and configuration. Standard deviations are calculated from saved repetition scalars with n−1 in the denominator. intNMF is a sensitivity experiment; the repaired TCGA run used separate missingness realizations."),
        ("Recovery failures", "Main benchmark recorded 38 ADNI and 9 TCGA refit failures. Imputation, refit, metric and baseline stages use saved status labels. Baseline failure slots for excluded TCGA CIA/O2PLS and original TCGA intNMF never reached refitting and are shown separately. Success rate uses attempted recovery runs; a metric-stage failure would be an unsuccessful run."),
        ("S7 pairing", "Both means use the same successful repetitions with finite latent RV and top-50 Jaccard. This yields 78 configurations per ADNI rate and 99 per TCGA rate. S7 covers the main benchmark."),
        ("Observed runtime", "Elapsed fitting plus evaluation excludes external imputation and screening. MI_SF timings use fully successful five-fit summaries divided by five. The timed denominator counts summaries; multiplying by fits per timed record gives the individual fits represented. Partial and failed records are excluded. Method-specific times are missing for TCGA MOFA and CMF."),
        ("Semi-synthetic runtime", "S9 pools mechanisms within dataset, method, removed modality and target rate. It uses positive elapsed seconds for successful refits and excludes shared imputation and subsequent metrics. intNMF sensitivity is outside S9."),
        ("Figure numbering", "Observed runtime appears in Figure 8. Semi-synthetic runtime appears in Figure S9."),
        ("Verification", f"{sum(checks.values()):,} numerical comparisons against saved records, figure tables and selected manuscript claims. See verification.txt for results and discrepancies."),
    ]
    if discrepancies:
        readme.extend(("Discrepancy", entry) for entry in discrepancies)
    tables = {"README": (["Topic", "Explanation"], readme), **tables}
    destination.mkdir(parents=True, exist_ok=True)
    output = workbook_bytes(tables)
    with targets[0].open("xb") as handle:
        handle.write(output)
    with targets[1].open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Worksheet", "Role", "CSV path relative to repository", "CSV SHA-256",
                         "Upstream path as recorded in repository provenance", "Upstream SHA-256"])
        writer.writerows(manifest)
    with targets[2].open("x", encoding="utf-8") as handle:
        handle.write("Final CSV and manuscript verification\n\n")
        for reference in (manuscript, supporting_pdf):
            if reference is not None:
                handle.write(f"{reference.name} SHA-256: {sha256(reference.read_bytes()).hexdigest()}\n")
        handle.write(f"{sum(checks.values()):,} numerical comparisons; {len(discrepancies)} numerical/provenance discrepancies.\n")
        for label, count in sorted(checks.items()):
            handle.write(f"{label}: {count}\n")
        for discrepancy in discrepancies:
            handle.write(f"DISCREPANCY: {discrepancy}\n")
        handle.write("\nWorkbook reopen and exact XML decimal checks passed for every worksheet.\n")
        for sheet, (_, rows) in tables.items():
            handle.write(f"{sheet}: {len(rows)} data rows\n")
    print(f"{targets[0]} | {sum(checks.values()):,} comparisons | {len(discrepancies)} discrepancies")


def workbook_bytes(tables):
    workbook = Workbook()
    workbook.remove(workbook.active)
    exact_cells = {}
    for name, (headers, rows) in tables.items():
        sheet = workbook.create_sheet(name)
        sheet.append(headers)
        exact = {}
        for row in rows:
            sheet.append([float(value) if isinstance(value, Decimal) else value for value in row])
            for column, value in enumerate(row, 1):
                if isinstance(value, Decimal):
                    cell = sheet.cell(sheet.max_row, column)
                    exact[cell.coordinate] = value
                    if headers[column - 1] in {"Target missingness rate", "Success rate"}:
                        cell.number_format = "0%" if headers[column - 1] == "Target missingness rate" else "0.00%"
            if sheet.max_row % 2 == 0:
                for cell in sheet[sheet.max_row]:
                    cell.fill = PatternFill("solid", fgColor="EDF5F8")
        sheet.freeze_panes = "C2" if name != "README" else "A2"
        sheet.auto_filter.ref = sheet.dimensions
        sheet.sheet_view.showGridLines = False
        sheet.row_dimensions[1].height = 32
        for column, header in enumerate(headers, 1):
            letter = get_column_letter(column)
            cell = sheet.cell(1, column)
            cell.fill = PatternFill("solid", fgColor="174760")
            cell.font = Font(color="FFFFFF", bold=True)
            cell.alignment = Alignment(vertical="center", wrap_text=True)
            width = max(len(header), *(len(str(row[column - 1])) for row in rows))
            sheet.column_dimensions[letter].width = min(59, max(12, width + 2))
        if name == "README":
            sheet.column_dimensions["B"].width = 104
            for i, row in enumerate(rows, 2):
                sheet.cell(i, 2).alignment = Alignment(wrap_text=True, vertical="center")
                sheet.row_dimensions[i].height = 17 * (1 + len(row[1]) // 95)
        exact_cells[f"xl/worksheets/sheet{len(exact_cells) + 1}.xml"] = exact

    # openpyxl's numeric writer shortens some decimals. Restore their source digits
    # while retaining the numeric cell type and checking every stored value.
    initial, result = BytesIO(), BytesIO()
    workbook.save(initial)
    pattern = re.compile(rb'(<c r="([A-Z]+[0-9]+)"[^>]* t="n"><v>)([^<]+)(</v>)')
    with ZipFile(initial) as source, ZipFile(result, "w") as export:
        for item in source.infolist():
            payload = source.read(item.filename)
            if item.filename in exact_cells:
                exact = exact_cells[item.filename]
                stored = set()

                def restore(match):
                    coordinate = match[2].decode()
                    if coordinate not in exact:
                        return match[0]
                    stored.add(coordinate)
                    return match[1] + str(exact[coordinate]).encode() + match[4]

                payload = pattern.sub(restore, payload)
                assert stored == exact.keys()
                for match in pattern.finditer(payload):
                    coordinate = match[2].decode()
                    if coordinate in exact:
                        assert Decimal(match[3].decode()).as_tuple() == exact[coordinate].as_tuple()
            export.writestr(item, payload)
    reopened = load_workbook(result, read_only=True)
    assert reopened.sheetnames == list(tables)
    for name, (headers, rows) in tables.items():
        actual_rows = list(reopened[name].values)
        assert actual_rows[0] == tuple(headers) and len(actual_rows) == len(rows) + 1
        for actual, expected in zip(actual_rows[1:], rows):
            assert actual == tuple(float(v) if isinstance(v, Decimal) else v for v in expected), (name, actual, expected)
    reopened.close()
    return result.getvalue()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "build/workbook")
    parser.add_argument("--manuscript", type=Path)
    parser.add_argument("--supporting-pdf", type=Path)
    args = parser.parse_args()
    main(args.output_dir, args.manuscript, args.supporting_pdf)
