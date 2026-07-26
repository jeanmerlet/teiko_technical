import sqlite3
import csv, os


ROOT = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(ROOT, "cell-count.csv")
SCHEMA_PATH = os.path.join(ROOT, "schema.sql")
DB_PATH = os.path.join(ROOT, "cell-count.db")


CELL_TYPES = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]


def main():
    # rebuild the database on each execution
    if os.path.exists(DB_PATH): os.remove(DB_PATH)

    # load csv
    with open(CSV_PATH, newline="") as csv_file:
        rows = list(csv.DictReader(csv_file))

    # csv dict to subjects set and samples and cell_counts lists
    subjects = {
        (
            row["subject"],
            row["project"],
            row["condition"],
            int(row["age"]),
            row["sex"],
            row["treatment"],
            row["response"] or None,
        )
        for row in rows
    }

    samples = [
        (
            row["sample"],
            row["subject"],
            row["sample_type"],
            int(row["time_from_treatment_start"]),
        )
        for row in rows
    ]

    cell_counts = [
        (row["sample"], cell_type, int(row[cell_type]))
        for row in rows
        for cell_type in CELL_TYPES
    ]

    # dict and lists to .db
    with sqlite3.connect(DB_PATH) as connection:
        with open(SCHEMA_PATH) as schema_file:
            connection.executescript(schema_file.read())

        connection.executemany(
            """
            INSERT INTO subjects
                (subject, project, condition, age, sex, treatment, response)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            subjects,
        )

        connection.executemany(
            """
            INSERT INTO samples
                (sample, subject, sample_type, time_from_treatment_start)
            VALUES (?, ?, ?, ?)
            """,
            samples,
        )

        connection.executemany(
            """
            INSERT INTO cell_counts (sample, cell_type, count)
            VALUES (?, ?, ?)
            """,
            cell_counts,
        )

    print(
        f"Created {os.path.basename(DB_PATH)} with "
        f"{len(subjects)} subjects, "
        f"{len(samples)} samples, and "
        f"{len(cell_counts)} cell counts."
    )


if __name__ == "__main__":
    main()
