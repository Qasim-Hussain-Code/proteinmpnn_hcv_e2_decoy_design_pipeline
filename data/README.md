# Data layout

`raw/` holds unchanged primary-source downloads and is ignored by Git. Retrieval manifests retain public URLs, accession snapshots and hashes. `processed/` holds the humanized experimental hybrid, mapping tables, discovery-only panels and the compact retained state models. `processed/quarantine/` contains held-out sequence data. Candidate-selection functions cannot consume it before a valid freeze. `work/` holds disposable modeling and design intermediates and is ignored.

Generated scientific outputs must never be replaced with synthetic fixtures. The offline fixture is explicitly synthetic and is used only for tests.
