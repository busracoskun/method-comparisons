"""Render manuscript figures from saved results and optional authorized UMAP coordinates."""
import argparse
from functools import partial
from pathlib import Path

import aggregate
import embeddings
from missingness_mechanisms import render_mechanisms
from observed_figures import clinical_data, clinical_figure, clinical_limits

ROOT = Path(__file__).resolve().parents[2]
OBSERVED = ROOT / "results/observed_data"
RECOVERY = ROOT / "results/semi_synthetic"
DIAGNOSTICS = ROOT / "results/diagnostics"

MAIN = {
    "1": "observed_modality_missingness",
    "2": "adni_clinical_alignment",
    "3": "tcga_clinical_alignment",
    "4": "observed_umap",
    "5": "adni_semi_synthetic_recovery",
    "6": "tcga_semi_synthetic_recovery",
    "7": "semi_synthetic_umap",
    "8": "observed_runtime",
}
SUPPLEMENTARY = {
    "s1": "adni_recovery_umap",
    "s2": "tcga_recovery_umap",
    "s3": "observed_umap",
    "s4": "feature_recovery",
    "s5": "adjusted_rand_index",
    "s6": "fitting_reliability",
    "s7": "patient_feature_recovery",
    "s8": "intnmf_sensitivity",
    "s9": "semi_synthetic_runtime",
    "s10": "missingness_mechanisms",
}
PRIVATE = {
    "4": embeddings.observed,
    "7": embeddings.recovery,
    "s1": partial(embeddings.recovery, dataset="ADNI"),
    "s2": partial(embeddings.recovery, dataset="TCGA"),
    "s3": partial(embeddings.observed, supplemental=True),
}
AGGREGATE = {
    "s4": (RECOVERY / "feature_recovery_plotted_values.csv", aggregate.feature),
    "s6": (DIAGNOSTICS / "fitting_reliability_plotted_values.csv", aggregate.reliability),
    "s7": (RECOVERY / "patient_feature_recovery_plotted_values.csv", aggregate.patient_feature),
    "s8": (RECOVERY / "intnmf_recovery_plotted_values.csv", aggregate.sensitivity),
    "s9": (DIAGNOSTICS / "semisynthetic_runtime_summary.csv", aggregate.runtime),
    "8": (DIAGNOSTICS / "observed_runtime_summary.csv", partial(aggregate.runtime, observed=True)),
}


def render(number, output, archive):
    if number in MAIN:
        target = output / "main" / f"figure_{int(number):02}_{MAIN[number]}"
    else:
        target = output / "supplementary" / f"supp_figure_{int(number[1:]):02}_{SUPPLEMENTARY[number]}"
    if number in PRIVATE:
        PRIVATE[number](target, archive)
    elif number in AGGREGATE:
        source, draw = AGGREGATE[number]
        draw(aggregate.read(source), target)
    elif number == "1":
        data = aggregate.read(OBSERVED / "missingness_by_group.csv")
        aggregate.missingness(data, target, hatch=False)
    elif number in {"2", "3"}:
        data, reference = clinical_data()
        cohort = "ADNI" if number == "2" else "TCGA"
        clinical_figure(data[data.dataset.eq(cohort)], reference, clinical_limits(data, reference), cohort, target)
    elif number in {"5", "6"}:
        data = aggregate.read(RECOVERY / "recovery_full_metric_summary.csv")
        aggregate.recovery(data, "ADNI" if number == "5" else "TCGA", target)
    elif number == "s5":
        data = aggregate.read(OBSERVED / "adjusted_rand_plotted_values.csv")
        reliability = aggregate.read(DIAGNOSTICS / "fitting_reliability_plotted_values.csv")
        data.dataset = data.dataset.str.upper()
        reliability.dataset = reliability.dataset.str.upper()
        aggregate.clinical(data, None, target, reliability=reliability.set_index(["dataset", "model", "style"]))
    elif number == "s10":
        render_mechanisms(RECOVERY, target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    all_figures = [*MAIN, *SUPPLEMENTARY]
    parser.add_argument("--figure", type=str.lower, choices=all_figures, action="append",
                        help="Select a figure by number. Repeat for several figures.")
    parser.add_argument("--output", type=Path, default=ROOT / "build/figures")
    parser.add_argument("--coordinate-archive", type=Path, help="Root of the authorized historical coordinate archive.")
    args = parser.parse_args()
    if args.figure and not args.coordinate_archive and any(number in PRIVATE for number in args.figure):
        parser.error("UMAP figures require --coordinate-archive.")
    output = args.output.resolve()
    if output.is_relative_to(ROOT / "figures") or output.is_relative_to(ROOT / "results"):
        parser.error("Choose an output directory outside the archived figures and results.")
    selected = args.figure or [n for n in all_figures if args.coordinate_archive or n not in PRIVATE]
    for number in selected:
        render(number, output, args.coordinate_archive)
    print(f"Rendered {len(selected)} figures in {output}")


if __name__ == "__main__":
    main()
