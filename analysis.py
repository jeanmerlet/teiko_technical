import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from statsmodels.formula.api import gee
from statsmodels.genmod.cov_struct import Exchangeable

import os


ROOT = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(ROOT, "cell-count.db")
OUTPUT_DIR = os.path.join(ROOT, "outputs")

FREQUENCIES_PATH = os.path.join(OUTPUT_DIR, "cell_frequencies.csv")
FILTERED_PATH = os.path.join(OUTPUT_DIR, "response_analysis_data.csv")
GEE_RESULTS_PATH = os.path.join(OUTPUT_DIR, "gee_results.csv")
SUBJECT_MEANS_PATH = os.path.join(OUTPUT_DIR, "subject_mean_frequencies.csv")
AGGREGATED_RESULTS_PATH = os.path.join(OUTPUT_DIR, "subject_aggregated_results.csv")
BOXPLOT_PATH = os.path.join(OUTPUT_DIR, "response_boxplots.png")

CELL_TYPES = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]


def load_database():
    with sqlite3.connect(DB_PATH) as connection:
        cell_counts = pd.read_sql_query(
            "SELECT sample, cell_type, count FROM cell_counts",
            connection,
        )

        metadata = pd.read_sql_query(
            """
            SELECT
                samples.sample,
                samples.subject,
                samples.sample_type,
                samples.time_from_treatment_start,
                subjects.condition,
                subjects.treatment,
                subjects.response
            FROM samples
            JOIN subjects
                ON samples.subject = subjects.subject
            """,
            connection,
        )

    return cell_counts, metadata


def calculate_frequencies(cell_counts):
    cell_counts["total_count"] = (
        cell_counts.groupby("sample")["count"].transform("sum")
    )

    cell_counts["percentage"] = (
        100 * cell_counts["count"] / cell_counts["total_count"]
    )

    frequencies = cell_counts.rename(
        columns={"cell_type": "population"}
    )[
        ["sample", "total_count", "population", "count", "percentage"]
    ]

    return frequencies


def get_response_data(frequencies, metadata):
    response_data = frequencies.merge(
        metadata,
        on="sample",
        validate="many_to_one",
    )

    response_data = response_data[
        (response_data["condition"] == "melanoma")
        & (response_data["treatment"] == "miraclib")
        & (response_data["sample_type"] == "PBMC")
        & (response_data["response"].isin(["yes", "no"]))
    ].copy()

    response_data["responder"] = (
        response_data["response"] == "yes"
    ).astype(int)

    return response_data


def run_gee_analysis(response_data):
    results = []

    for population in CELL_TYPES:
        population_data = response_data[
            response_data["population"] == population
        ].copy()

        model = gee(
            "percentage ~ responder + C(time_from_treatment_start)",
            groups="subject",
            data=population_data,
            family=sm.families.Gaussian(),
            cov_struct=Exchangeable(),
        )

        fit = model.fit()

        responders = population_data[
            population_data["response"] == "yes"
        ]

        nonresponders = population_data[
            population_data["response"] == "no"
        ]

        confidence_interval = fit.conf_int().loc["responder"]

        results.append(
            {
                "population": population,
                "responder_subjects": responders["subject"].nunique(),
                "nonresponder_subjects": (
                    nonresponders["subject"].nunique()
                ),
                "responder_mean_percentage": (
                    responders["percentage"].mean()
                ),
                "nonresponder_mean_percentage": (
                    nonresponders["percentage"].mean()
                ),
                "effect": fit.params["responder"],
                "standard_error": fit.bse["responder"],
                "confidence_interval_lower": confidence_interval.iloc[0],
                "confidence_interval_upper": confidence_interval.iloc[1],
                "p_value": fit.pvalues["responder"],
            }
        )

    results = pd.DataFrame(results)

    results["adjusted_p_value"] = multipletests(
        results["p_value"],
        method="fdr_bh",
    )[1]

    results["significant"] = results["adjusted_p_value"] < 0.05

    return results


def run_subject_aggregated_analysis(response_data):
    subject_means = (
        response_data.groupby(
            ["subject", "response", "population"],
            as_index=False,
        )["percentage"]
        .mean()
    )

    results = []

    for population in CELL_TYPES:
        population_data = subject_means[
            subject_means["population"] == population
        ]

        responders = population_data.loc[
            population_data["response"] == "yes",
            "percentage",
        ]

        nonresponders = population_data.loc[
            population_data["response"] == "no",
            "percentage",
        ]

        test = mannwhitneyu(
            responders,
            nonresponders,
            alternative="two-sided",
        )

        results.append(
            {
                "population": population,
                "responder_subjects": len(responders),
                "nonresponder_subjects": len(nonresponders),
                "responder_median_percentage": responders.median(),
                "nonresponder_median_percentage": nonresponders.median(),
                "median_difference": (
                    responders.median() - nonresponders.median()
                ),
                "test_statistic": test.statistic,
                "p_value": test.pvalue,
            }
        )

    results = pd.DataFrame(results)

    results["adjusted_p_value"] = multipletests(
        results["p_value"],
        method="fdr_bh",
    )[1]

    results["significant"] = results["adjusted_p_value"] < 0.05

    return subject_means, results


def plot_response_groups(response_data):
    plt.figure(figsize=(11, 6))

    sns.boxplot(
        data=response_data,
        x="population",
        y="percentage",
        hue="response",
        order=CELL_TYPES,
        hue_order=["no", "yes"],
    )

    plt.xlabel("Immune cell population")
    plt.ylabel("Relative frequency (%)")
    plt.title(
        "Immune Cell Frequencies in Miraclib-Treated Melanoma Patients"
    )
    plt.legend(title="Response")
    plt.tight_layout()
    plt.savefig(BOXPLOT_PATH, dpi=300)
    plt.close()


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    cell_counts, metadata = load_database()

    frequencies = calculate_frequencies(cell_counts)
    frequencies.to_csv(FREQUENCIES_PATH, index=False)
    response_data = get_response_data(frequencies, metadata)
    response_data.to_csv(FILTERED_PATH, index=False)
    gee_results = run_gee_analysis(response_data)
    gee_results.to_csv(GEE_RESULTS_PATH, index=False)
    subject_means, aggregated_results = run_subject_aggregated_analysis(response_data)
    subject_means.to_csv(SUBJECT_MEANS_PATH, index=False)
    aggregated_results.to_csv(AGGREGATED_RESULTS_PATH, index=False)
    plot_response_groups(response_data)

    print("\nGEE results:")
    print(gee_results.to_string(index=False))
    print("\nSubject-aggregated results:")
    print(aggregated_results.to_string(index=False))
    print(f"\nSaved analysis outputs to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
