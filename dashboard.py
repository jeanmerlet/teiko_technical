import pandas as pd
import plotly.express as px
import streamlit as st

import os


ROOT = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(ROOT, "outputs")

FREQUENCIES_PATH = os.path.join(OUTPUT_DIR, "cell_frequencies.csv")
RESPONSE_DATA_PATH = os.path.join(OUTPUT_DIR, "response_analysis_data.csv")
GEE_RESULTS_PATH = os.path.join(OUTPUT_DIR, "gee_results.csv")
AGGREGATED_RESULTS_PATH = os.path.join(OUTPUT_DIR, "subject_aggregated_results.csv")
BASELINE_SAMPLES_PATH = os.path.join(OUTPUT_DIR, "baseline_melanoma_pbmc_samples.csv")
BASELINE_PROJECTS_PATH = os.path.join(OUTPUT_DIR, "baseline_samples_by_project.csv")
BASELINE_RESPONSE_PATH = os.path.join(OUTPUT_DIR, "baseline_subjects_by_response.csv")
BASELINE_SEX_PATH = os.path.join(OUTPUT_DIR, "baseline_subjects_by_sex.csv")
B_CELL_AVERAGE_PATH = os.path.join(OUTPUT_DIR, "melanoma_male_responder_b_cell_average.csv")

CELL_TYPES = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]


st.set_page_config(
    page_title="Drug Candidate Analysis",
    page_icon="🧬",
    layout="wide",
)


@st.cache_data
def load_outputs():
    paths = [
        FREQUENCIES_PATH,
        RESPONSE_DATA_PATH,
        GEE_RESULTS_PATH,
        AGGREGATED_RESULTS_PATH,
        BASELINE_SAMPLES_PATH,
        BASELINE_PROJECTS_PATH,
        BASELINE_RESPONSE_PATH,
        BASELINE_SEX_PATH,
        B_CELL_AVERAGE_PATH,
    ]

    missing = [path for path in paths if not os.path.exists(path)]

    if missing:
        st.error("Required output files are missing. Run `make pipeline` before starting the dashboard.")
        st.stop()

    return {
        "frequencies": pd.read_csv(FREQUENCIES_PATH),
        "response_data": pd.read_csv(RESPONSE_DATA_PATH),
        "gee_results": pd.read_csv(GEE_RESULTS_PATH),
        "aggregated_results": pd.read_csv(AGGREGATED_RESULTS_PATH),
        "baseline_samples": pd.read_csv(BASELINE_SAMPLES_PATH),
        "samples_by_project": pd.read_csv(BASELINE_PROJECTS_PATH),
        "subjects_by_response": pd.read_csv(BASELINE_RESPONSE_PATH),
        "subjects_by_sex": pd.read_csv(BASELINE_SEX_PATH),
    }


outputs = load_outputs()

frequencies = outputs["frequencies"]
response_data = outputs["response_data"]
gee_results = outputs["gee_results"]
aggregated_results = outputs["aggregated_results"]
baseline_samples = outputs["baseline_samples"]
samples_by_project = outputs["samples_by_project"]
subjects_by_response = outputs["subjects_by_response"]
subjects_by_sex = outputs["subjects_by_sex"]


st.title("Drug Candidate Analysis")

overview_tab, response_tab, baseline_tab = st.tabs(
    [
        "Data Overview",
        "Miraclib Response",
        "Baseline Samples",
    ]
)


with overview_tab:
    st.header("Cell population frequencies")
    st.write(
        "The table on this tab shows per-sample cell type frequencies for all samples, with these columns:"
    )
    st.markdown("""
        - sample: sample id
        - total_count: total cell count
        - population: cell type
        - count: population cell count
        - percentage: population relative frequency
    """)

    selected_populations = st.multiselect(
        "Cell types",
        CELL_TYPES,
        default=CELL_TYPES
    )

    sample_search = st.text_input(
        "Sample",
        placeholder="Enter all or part of a sample ID"
    )

    frequency_display = frequencies[
        frequencies["population"].isin(selected_populations)
    ].copy()

    if sample_search:
        frequency_display = frequency_display[
            frequency_display["sample"].str.contains(
                sample_search,
                case=False,
                regex=False,
            )
        ]

    frequency_display["percentage"] = frequency_display["percentage"].round(2)

    st.dataframe(
        frequency_display,
        use_container_width=True,
        hide_index=True,
    )

    st.download_button(
        "Download filtered frequencies",
        frequency_display.to_csv(index=False),
        file_name="filtered_cell_frequencies.csv",
        mime="text/csv",
    )


with response_tab:
    st.header("Miraclib response analysis")

    st.write(
        "The box plots below are generated with the data only from the miraclib-treated melanoma PBMC samples. Generalized estimating equations (GEE) identified a small but statistically significant 0.64-percentage-point increase in CD4+ T-cell frequency among responders compared with nonresponders (30.54% versus 29.90%). No other cell populations differed significantly after multiple-testing correction. Subject-aggregated Welch’s t-tests produced consistent results, whereas the nonparametric Mann–Whitney U sensitivity analysis identified no significant differences after correction."
    )

    available_times = sorted(response_data["time_from_treatment_start"].unique())

    selected_times = st.multiselect(
        "Treatment times",
        available_times,
        default=available_times,
        format_func=lambda time: f"Day {time}",
    )

    selected_response_populations = st.multiselect(
        "Populations",
        CELL_TYPES,
        default=CELL_TYPES,
        key="response_populations",
    )

    plot_data = response_data[
        response_data["time_from_treatment_start"].isin(selected_times)
        & response_data["population"].isin(selected_response_populations)
    ]

    figure = px.box(
        plot_data,
        x="population",
        y="percentage",
        color="response",
        category_orders={
            "population": CELL_TYPES,
            "response": ["no", "yes"],
        },
        color_discrete_map={
            "no": "#377eb8",
            "yes": "#e67e22",
        },
        labels={
            "population": "Immune cell population",
            "percentage": "Relative frequency (%)",
            "response": "Response",
        },
        title="Relative Cell Frequencies by Treatment Response",
    )

    st.plotly_chart(figure, use_container_width=True)

    st.subheader("GEE results")


    gee_display = gee_results[
        [
            "population",
            "effect",
            "confidence_interval_lower",
            "confidence_interval_upper",
            "adjusted_p_value",
            "significant",
        ]
    ].copy()

    gee_display.columns = [
        "Population",
        "Effect",
        "95% CI lower",
        "95% CI upper",
        "FDR-adjusted p-value",
        "Significant",
    ]

    gee_display[["Effect", "95% CI lower", "95% CI upper"]] = gee_display[
        ["Effect", "95% CI lower", "95% CI upper"]
    ].round(3)

    gee_display["FDR-adjusted p-value"] = gee_display[
        "FDR-adjusted p-value"
    ].round(4)

    st.dataframe(
        gee_display,
        use_container_width=True,
        hide_index=True,
    )

    with st.expander("Subject-aggregated sensitivity analysis"):
        aggregated_display = aggregated_results[
            [
                "population",
                "mean_difference",
                "welch_adjusted_p_value",
                "welch_significant",
                "mann_whitney_adjusted_p_value",
                "mann_whitney_significant",
            ]
        ].copy()

        aggregated_display.columns = [
            "Population",
            "Mean difference",
            "Welch adjusted p-value",
            "Welch significant",
            "Mann–Whitney adjusted p-value",
            "Mann–Whitney significant",
        ]

        aggregated_display = aggregated_display.round(4)

        st.dataframe(
            aggregated_display,
            use_container_width=True,
            hide_index=True,
        )


with baseline_tab:
    st.header("Miraclib-treated subset")

    st.write(
        "The table on this tab is filtered to baseline PBMC samples from melanoma patients treated with miraclib."
    )

    project_column, response_column, sex_column = st.columns(3)

    with project_column:
        st.subheader("Samples by project")
        st.dataframe(
            samples_by_project,
            use_container_width=True,
            hide_index=True,
        )

    with response_column:
        st.subheader("Subjects by response")
        st.dataframe(
            subjects_by_response,
            use_container_width=True,
            hide_index=True,
        )

    with sex_column:
        st.subheader("Subjects by sex")
        st.dataframe(
            subjects_by_sex,
            use_container_width=True,
            hide_index=True,
        )

    st.subheader("Qualifying baseline samples")

    st.dataframe(
        baseline_samples,
        use_container_width=True,
        hide_index=True,
    )

    st.download_button(
        "Download baseline samples",
        baseline_samples.to_csv(index=False),
        file_name="baseline_melanoma_pbmc_samples.csv",
        mime="text/csv",
    )
