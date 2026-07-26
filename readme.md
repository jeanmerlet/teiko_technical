# Teiko Technical Assessment

This repository contains a reproducible SQLite and Python pipeline for loading, analyzing, and visualizing the provided immune cell-count dataset.

Dashboard: **[link to be added after deployment]**

## Running the project

The project is designed to run in GitHub Codespaces.

Install the required Python dependencies:

```bash
make setup
```

Initialize the SQLite database and generate all analysis outputs:

```bash
make pipeline
```

Start the interactive dashboard:

```bash
make dashboard
```

The pipeline can also be run directly:

```bash
python load_data.py
python analysis.py
```

## Database schema

The SQLite database uses three related tables:

* `subjects` contains subject-level metadata, including project, condition, age, sex, treatment, and response.
* `samples` contains one row per biological sample, with its subject, sample type, and collection time.
* `cell_counts` contains one row per sample and immune cell population.

Each sample references a subject through a foreign key, and each cell-count record references a sample. The composite primary key `(sample, cell_type)` ensures that each sample has at most one count for each population.

This design avoids repeating subject and sample metadata for every cell population. Storing cell counts in long format also allows additional populations (e.g., granulocytes) to be added without changing the database schema.

For hundreds of projects and thousands of samples, the same schema would remain suitable. Indexes could be added to frequently filtered columns such as project, condition, treatment, and collection time. Larger aggregations could be performed directly in SQL so that only summarized results are loaded into Python. If data volume or concurrent access grew beyond SQLite’s intended use, the schema could be transferred to a server-based relational database such as PostgreSQL.

## Code structure

* `cell-count.csv`: provided input data.
* `schema.sql`: SQLite table definitions, relationships, and data constraints.
* `load_data.py`: rebuilds the database and loads all CSV records.
* `analysis.py`: calculates cell frequencies, performs response-group statistical analyses, creates plots, and runs the requested baseline subset queries.
* `dashboard.py`: interactive display of the generated tables and plots.
* `outputs/`: generated analysis tables and figures.
* `Makefile`: commands for dependency installation, pipeline execution, and dashboard startup.

The loading, analysis, and dashboard components are kept separate so that the database and all outputs can be reproduced without starting the dashboard.

## Analysis

Relative cell frequencies are calculated by dividing each population count by the total count for its sample.

Miraclib-treated melanoma PBMC samples are compared between responders and nonresponders using generalized estimating equations. The models adjust for collection time and account for repeated samples from the same subject. Subject-aggregated Welch’s t-tests are included as a sensitivity analysis, with Mann–Whitney U tests as a nonparametric check. Benjamini–Hochberg correction is applied across the five immune populations.

The pipeline also identifies baseline melanoma PBMC samples from miraclib-treated subjects and reports sample counts by project and subject counts by response and sex.
