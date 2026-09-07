# Result layout

Both operations return a deep replacement of the APB2 result. Neither mutates its input, and both refuse to run when an output name already exists.

## Peptide verification

Writes a feature-aligned `varm["fasta_validation"]` table on every peptide-derived level it covered.

| Column | Type | Meaning |
| --- | --- | --- |
| `peptide_in_fasta` | Boolean | The stripped sequence occurs in at least one protein |
| `fasta_match_site_count` | UInt64 | Total match sites across the database |
| `fasta_matching_protein_count` | UInt64 | Distinct proteins containing the sequence |
| `fasta_matching_protein_ids` | String | Separator-joined matching accessions |
| `reported_member_count` | UInt64 | Protein-group members the result reported |
| `reported_members_in_fasta_count` | UInt64 | Reported members found in the database |
| `peptide_in_reported_protein` | Boolean | The sequence occurs in a reported member |

The last three columns are null when the level reports no protein group, which distinguishes "not reported" from "reported and unmatched".

## Protein annotation

Writes three things:

- `varm["fasta"]` on the `protein` level — one row per protein group
- annotation table `fasta_protein_group_members` — one row per group member, lossless
- feature relation `fasta_protein_group_membership` — directed `member_of` mapping from that table to the `protein` axis

`varm["fasta"]` columns:

| Column | Type |
| --- | --- |
| `reported_member_count`, `matched_member_count`, `unmatched_member_count`, `ambiguous_member_count` | UInt64 |
| `all_members_in_fasta`, `any_member_in_fasta` | Boolean |
| `fasta_descriptions`, `fasta_gene_names`, `fasta_organism_names` | String |

The member table keeps `protein_group`, `protein_member`, `member_ordinal`, `match_ordinal`, `match_status`, `fasta_id`, `fasta_description`, and `fasta_sequence_length`. A group is never collapsed to its leading accession — the group stays one row on the protein axis while every member keeps its own row in the annotation table.

## Provenance

Each operation records itself under result metadata key `fasta`, as `peptide_verification` or `protein_annotation`, together with the FASTA source paths, checksums and ordinals, the separator, and — for verification — the requested and resolved Prozor backend. Applying the same operation twice is refused rather than silently re-recorded.

## Required protein-frame columns

`protein_fasta` must supply `id`, `description`, `sequence`, `fasta_source_path`, `fasta_source_checksum`, `fasta_source_ordinal`, and `fasta_record_ordinal`. A null `id` or `sequence` is an error.
