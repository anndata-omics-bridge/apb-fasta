# Result layout

Every operation returns a deep replacement of the APB2 result. None mutates its input, and each refuses to run when an output name already exists.

## Peptide verification

Writes a feature-aligned `varm["fasta_validation"]` table on every peptide-derived level it covered.

| Column | Type | Meaning |
| --- | --- | --- |
| `peptide_in_fasta` | Boolean | The stripped sequence occurs in at least one protein |
| `fasta_match_site_count` | UInt64 | Total match sites across the database |
| `fasta_matching_protein_count` | UInt64 | Distinct proteins containing the sequence |
| `fasta_matching_protein_ids` | String | `;`-joined FASTA ids of the matching proteins |
| `fasta_matching_organisms` | String | Sorted, `;`-joined organism mnemonics of the matching proteins, e.g. `HUMAN;YEAST` |
| `fasta_matches_contaminant` | Boolean | A matching protein is a contaminant (`is_contaminant`) |
| `reported_member_count` | UInt64 | Protein-group members the result reported |
| `reported_members_in_fasta_count` | UInt64 | Reported members found in the database |
| `reported_member_organisms` | String | Sorted, `;`-joined organism mnemonics of the reported members found |
| `reported_members_contaminant` | Boolean | A reported member found is a contaminant |
| `peptide_in_reported_protein` | Boolean | The sequence occurs in a reported member |
| `fasta_il_only` | Boolean | `peptide_in_fasta` or `peptide_in_reported_protein` holds only for another I/L spelling |
| `is_protein_group_specific` | Boolean | Every protein the peptide matches is a member of its own protein group (`protein_assignment`) |
| `is_proteotypic` | Boolean | The peptide matches exactly one protein of the FASTA |
| `reported_leading_id` | String | FASTA `id` of the leading protein, e.g. `sp\|P02768\|ALBU_HUMAN` |
| `reported_leading_accession` | String | Its accession |
| `reported_leading_description` | String | Its header text after the id |
| `reported_leading_gene_name` | String | Its gene name |
| `reported_leading_protein_length` | UInt64 | Its sequence length |
| `reported_leading_tryptic_peptides` | UInt64 | Its tryptic peptides of 7 to 29 residues, cleaved after K or R unless P follows |

The four `reported_member*` columns and `peptide_in_reported_protein` are null when the level reports no protein group, which distinguishes "not reported" from "reported and unmatched". `is_protein_group_specific` is null when the peptide matches no protein, the level declares no `protein_assignment`, or no member of the feature's group is in the FASTA. It tests the `protein_assignment` group, the vendor's inferred group, not the `fasta_accessions` list the `reported_*` columns read, which can name every protein a peptide maps to. `is_proteotypic` is null only when the peptide matches no protein. The two differ: a peptide shared by two proteins of one group is group specific but not proteotypic.

The leading protein is the first member of the level's `protein_assignment` column, found as any reported member is; the first such record in the FASTA wins. Its six `reported_leading_*` columns carry what prolfquapp reads from a FASTA for a protein, counted as prolfquapp's `get_annot_from_fasta` counts; they are null when the level declares no `protein_assignment` or the FASTA lacks that member.

Reported members are read from apb2's `;`-separated protein lists. A member names a FASTA record by its id, accession or entry name, or, for a `db|ACC|NAME` id, by `ACC|NAME` as PEAKS writes it or by `ACC`, which keeps a contaminant's `Cont_` prefix.

Mass spectrometry cannot tell isoleucine from leucine, so a vendor may spell a peptide differently from its FASTA protein. Protein sequences are matched as read. A peptide that occurs exactly keeps exactly its proteins; its other I/L spellings are tried only when it occurs nowhere, or in no reported member the database contains. In the first case the spellings' matches become its proteins; in the second they count only toward `peptide_in_reported_protein`. Coverage counts the targets matched through another spelling as `il_only_matched_feature_count`.

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

## Peptide properties

`add_peptide_properties(parsed)` needs no FASTA. It writes a feature-aligned `varm["peptide_properties"]` table on every peptide-derived level, computed from `ProForma_peptide` by protein_fasta's `peptide_property_frame()`: length, average molecular weight, isoelectric point (EMBOSS), Kyte-Doolittle hydrophobicity, instability and Boman indices, charge at pH 7 (Sillero), predicted reversed-phase retention time (Goloborodko et al. 2010 coefficients), missed cleavages, proline count, C-terminal residue, and Cys/Met/Trp, N-terminal Q/E or Cys, NG and DP flags. The protein_fasta API reference lists the columns and their types. A feature without a sequence, or whose sequence has a residue outside the 20 standard amino acids, has null properties. The installed protein_fasta version is recorded under `fasta.provenance.peptide_properties`.

## Provenance

FASTA metadata schema 5 groups root sources and settings by operation. `fasta.provenance.peptide_verification` and `fasta.provenance.protein_annotation` hold source file names, checksums, ordinals and database identity; peptide verification also records the requested/resolved Prozor backend. `fasta.provenance.peptide_properties` holds only the installed `protein_fasta_version`. Each checked level carries its own `fasta` record: `result.peptide_verification` on peptide-derived levels and `result.protein_annotation` on the protein level, a `summary` of matched, unmatched and decoy peptides or matched, unmatched and ambiguous members, and `details` naming the `varm` tables, member table and relation that hold the per-feature evidence. MuData stores provenance once in its root `uns["apb"]["fasta"]` and each level's record in that modality; a standalone H5AD keeps the same two parts, with the level's under `uns[<level>]["apb"]`. Applying the same operation twice is refused rather than silently re-recorded.

## Required protein-frame columns

`protein_fasta` must supply `id`, `description`, `sequence`, `fasta_source_path`, `fasta_source_checksum`, `fasta_source_ordinal`, and `fasta_record_ordinal`. A null `id` or `sequence` is an error.
