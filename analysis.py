import sqlite3
import pandas as pd
import os


ROOT = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(ROOT, "cell-count.db")
OUTPUT_DIR = os.path.join(ROOT, "outputs")
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "cell_frequencies.csv")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with sqlite3.connect(DB_PATH) as connection:
        cell_counts = pd.read_sql_query(
            "SELECT sample, cell_type, count FROM cell_counts",
            connection,
        )

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

    frequencies.to_csv(OUTPUT_PATH, index=False)

    print(frequencies.head(10).to_string(index=False))
    print(f"\nSaved {len(frequencies)} rows to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
