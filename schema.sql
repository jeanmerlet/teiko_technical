PRAGMA foreign_keys = ON;

CREATE TABLE subjects (
    subject TEXT PRIMARY KEY,
    project TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER NOT NULL CHECK (age >= 0),
    sex TEXT NOT NULL,
    treatment TEXT NOT NULL,
    response TEXT CHECK (response IN ('yes', 'no') OR response IS NULL)
);

CREATE TABLE samples (
    sample TEXT PRIMARY KEY,
    subject TEXT NOT NULL,
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    FOREIGN KEY (subject) REFERENCES subjects(subject)
);

CREATE TABLE cell_counts (
    sample TEXT NOT NULL,
    cell_type TEXT NOT NULL CHECK (
        cell_type IN (
            'b_cell',
            'cd8_t_cell',
            'cd4_t_cell',
            'nk_cell',
            'monocyte'
        )
    ),
    count INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample, cell_type),
    FOREIGN KEY (sample) REFERENCES samples(sample)
);
