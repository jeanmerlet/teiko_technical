import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm

from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.formula.api import gee
from statsmodels.genmod.cov_struct import Exchangeable
from statsmodels.stats.multitest import multipletests

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
TIME_BOXPLOT_PATH = os.path.join(OUTPUT_DIR, "response_boxplots_by_time.png")
BASELINE_SAMPLES_PATH = os.path.join(OUTPUT_DIR, "baseline_melanoma_pbmc_samples.csv")
BASELINE_PROJECTS_PATH = os.path.join(OUTPUT_DIR, "baseline_samples_by_project.csv")
BASELINE_RESPONSE_PATH = os.path.join(OUTPUT_DIR, "baseline_subjects_by_response.csv")
BASELINE_SEX_PATH = os.path.join(OUTPUT_DIR, "baseline_subjects_by_sex.csv")
MELANOMA_MALE_B_CELL_PATH = os.path.join(OUTPUT_DIR, "melanoma_male_responder_b_cell_average.csv")

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
    cell_counts["total_count"] = cell_counts.groupby("sample")["count"].transform("sum")
    cell_counts["percentage"] = 100 * cell_counts["count"] / cell_counts["total_count"]

    frequencies = cell_counts.rename(columns={"cell_type": "population"})[
        ["sample", "total_count", "population", "count", "percentage"]
    ]

    return frequencies


def get_response_data(frequencies, metadata):
    response_data = frequencies.merge(metadata, on="sample", validate="many_to_one")

    response_data = response_data[
        (response_data["condition"] == "melanoma")
        & (response_data["treatment"] == "miraclib")
        & (response_data["sample_type"] == "PBMC")
        & (response_data["response"].isin(["yes", "no"]))
    ].copy()

    response_data["responder"] = (response_data["response"] == "yes").astype(int)

    return response_data


def run_gee_analysis(response_data):
    results = []

    for population in CELL_TYPES:
        population_data = response_data[response_data["population"] == population].copy()

        model = gee(
            "percentage ~ responder + C(time_from_treatment_start)",
            groups="subject",
            data=population_data,
            family=sm.families.Gaussian(),
            cov_struct=Exchangeable(),
        )

        fit = model.fit()
        responders = population_data[population_data["response"] == "yes"]
        nonresponders = population_data[population_data["response"] == "no"]
        confidence_interval = fit.conf_int().loc["responder"]

        results.append(
            {
                "population": population,
                "responder_subjects": responders["subject"].nunique(),
                "nonresponder_subjects": nonresponders["subject"].nunique(),
                "responder_mean_percentage": responders["percentage"].mean(),
                "nonresponder_mean_percentage": nonresponders["percentage"].mean(),
                "effect": fit.params["responder"],
                "standard_error": fit.bse["responder"],
                "confidence_interval_lower": confidence_interval.iloc[0],
                "confidence_interval_upper": confidence_interval.iloc[1],
                "p_value": fit.pvalues["responder"],
            }
        )

    results = pd.DataFrame(results)
    results["adjusted_p_value"] = multipletests(results["p_value"], method="fdr_bh")[1]
    results["significant"] = results["adjusted_p_value"] < 0.05

    return results


def run_subject_aggregated_analysis(response_data):
    subject_means = response_data.groupby(
        ["subject", "response", "population"],
        as_index=False,
    )["percentage"].mean()

    results = []

    for population in CELL_TYPES:
        population_data = subject_means[subject_means["population"] == population]
        responders = population_data.loc[population_data["response"] == "yes", "percentage"]
        nonresponders = population_data.loc[population_data["response"] == "no", "percentage"]

        welch_test = ttest_ind(responders, nonresponders, equal_var=False)
        mann_whitney_test = mannwhitneyu(responders, nonresponders, alternative="two-sided")

        results.append(
            {
                "population": population,
                "responder_subjects": len(responders),
                "nonresponder_subjects": len(nonresponders),
                "responder_mean_percentage": responders.mean(),
                "nonresponder_mean_percentage": nonresponders.mean(),
                "mean_difference": responders.mean() - nonresponders.mean(),
                "welch_t_statistic": welch_test.statistic,
                "welch_p_value": welch_test.pvalue,
                "responder_median_percentage": responders.median(),
                "nonresponder_median_percentage": nonresponders.median(),
                "median_difference": responders.median() - nonresponders.median(),
                "mann_whitney_statistic": mann_whitney_test.statistic,
                "mann_whitney_p_value": mann_whitney_test.pvalue,
            }
        )

    results = pd.DataFrame(results)
    results["welch_adjusted_p_value"] = multipletests(results["welch_p_value"], method="fdr_bh")[1]
    results["welch_significant"] = results["welch_adjusted_p_value"] < 0.05
    results["mann_whitney_adjusted_p_value"] = multipletests(results["mann_whitney_p_value"], method="fdr_bh")[1]
    results["mann_whitney_significant"] = results["mann_whitney_adjusted_p_value"] < 0.05

    return subject_means, results


def run_baseline_subset_analysis():
    with sqlite3.connect(DB_PATH) as connection:
        baseline_samples = pd.read_sql_query(
            """
            SELECT
                samples.sample,
                samples.subject,
                subjects.project,
                subjects.response,
                subjects.sex
            FROM samples
            JOIN subjects
                ON samples.subject = subjects.subject
            WHERE subjects.condition = 'melanoma'
                AND samples.sample_type = 'PBMC'
                AND subjects.treatment = 'miraclib'
                AND samples.time_from_treatment_start = 0
            ORDER BY samples.sample
            """,
            connection,
        )

        samples_by_project = pd.read_sql_query(
            """
            SELECT
                subjects.project,
                COUNT(samples.sample) AS sample_count
            FROM samples
            JOIN subjects
                ON samples.subject = subjects.subject
            WHERE subjects.condition = 'melanoma'
                AND samples.sample_type = 'PBMC'
                AND subjects.treatment = 'miraclib'
                AND samples.time_from_treatment_start = 0
            GROUP BY subjects.project
            ORDER BY subjects.project
            """,
            connection,
        )

        subjects_by_response = pd.read_sql_query(
            """
            SELECT
                subjects.response,
                COUNT(DISTINCT samples.subject) AS subject_count
            FROM samples
            JOIN subjects
                ON samples.subject = subjects.subject
            WHERE subjects.condition = 'melanoma'
                AND samples.sample_type = 'PBMC'
                AND subjects.treatment = 'miraclib'
                AND samples.time_from_treatment_start = 0
                AND subjects.response IN ('yes', 'no')
            GROUP BY subjects.response
            ORDER BY subjects.response
            """,
            connection,
        )

        subjects_by_sex = pd.read_sql_query(
            """
            SELECT
                CASE subjects.sex
                    WHEN 'M' THEN 'male'
                    WHEN 'F' THEN 'female'
                    ELSE subjects.sex
                END AS sex,
                COUNT(DISTINCT samples.subject) AS subject_count
            FROM samples
            JOIN subjects
                ON samples.subject = subjects.subject
            WHERE subjects.condition = 'melanoma'
                AND samples.sample_type = 'PBMC'
                AND subjects.treatment = 'miraclib'
                AND samples.time_from_treatment_start = 0
            GROUP BY subjects.sex
            ORDER BY subjects.sex
            """,
            connection,
        )

    return baseline_samples, samples_by_project, subjects_by_response, subjects_by_sex


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
    plt.title("Immune Cell Frequencies in Miraclib-Treated Melanoma Patients")
    plt.legend(title="Response")
    plt.tight_layout()
    plt.savefig(BOXPLOT_PATH, dpi=300)
    plt.close()


def plot_response_groups_by_time(response_data):
    plot = sns.catplot(
        data=response_data,
        x="population",
        y="percentage",
        hue="response",
        col="time_from_treatment_start",
        kind="box",
        order=CELL_TYPES,
        hue_order=["no", "yes"],
        height=5,
        aspect=1,
        sharey=True,
    )

    plot.set_axis_labels("Immune cell population", "Relative frequency (%)")
    plot.set_titles("Day {col_name}")
    plot.fig.suptitle("Immune Cell Frequencies by Treatment Time", y=1.04)
    plot.savefig(TIME_BOXPLOT_PATH, dpi=300, bbox_inches="tight")
    plt.close(plot.fig)


def get_melanoma_male_b_cell_average():
    with sqlite3.connect(DB_PATH) as connection:
        result = pd.read_sql_query(
            """
            SELECT
                ROUND(AVG(cell_counts.count), 2) AS average_b_cell_count
            FROM cell_counts
            JOIN samples
                ON cell_counts.sample = samples.sample
            JOIN subjects
                ON samples.subject = subjects.subject
            WHERE subjects.condition = 'melanoma'
                AND subjects.sex = 'M'
                AND subjects.response = 'yes'
                AND samples.time_from_treatment_start = 0
                AND cell_counts.cell_type = 'b_cell'
            """,
            connection,
        )

    return result


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

    baseline_samples, samples_by_project, subjects_by_response, subjects_by_sex = run_baseline_subset_analysis()
    baseline_samples.to_csv(BASELINE_SAMPLES_PATH, index=False)
    samples_by_project.to_csv(BASELINE_PROJECTS_PATH, index=False)
    subjects_by_response.to_csv(BASELINE_RESPONSE_PATH, index=False)
    subjects_by_sex.to_csv(BASELINE_SEX_PATH, index=False)

    plot_response_groups(response_data)
    plot_response_groups_by_time(response_data)

    b_cell_average = get_melanoma_male_b_cell_average()
    b_cell_average.to_csv(MELANOMA_MALE_B_CELL_PATH, index=False)

    gee_print = gee_results[["population", "effect", "adjusted_p_value", "significant"]]
    aggregated_print = aggregated_results[["population", "mean_difference", "welch_adjusted_p_value", "welch_significant"]]
    average = b_cell_average.loc[0, "average_b_cell_count"]

    print("\nGEE results:\n", gee_print.to_string(index=False))
    print("\nSubject-aggregated Welch results:\n", aggregated_print.to_string(index=False))
    print("\nBaseline samples by project:\n", samples_by_project.to_string(index=False))
    print("\nBaseline subjects by response:\n", subjects_by_response.to_string(index=False))
    print("\nBaseline subjects by sex:\n", subjects_by_sex.to_string(index=False))
    print(f"\nSaved analysis outputs to {OUTPUT_DIR}")
    print(f"\nAverage B-cell count: {average:.2f}")


if __name__ == "__main__":
    main()
