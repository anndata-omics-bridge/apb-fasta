"""User-selected FASTA annotation behavior."""

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True, slots=True)
class FastaAnnotationParameters:
    """Configuration shared by peptide validation and protein annotation."""

    protein_group_separator: str = ";"
    matcher_backend: Literal["auto", "ahocorapy", "ahocorasick_rs"] = "auto"
    il_equivalent: bool = False

    def __post_init__(self) -> None:
        """Reject a separator that cannot define protein-group members."""
        if not self.protein_group_separator:
            raise ValueError("protein_group_separator must not be empty")


DEFAULT_FASTA_ANNOTATION_PARAMETERS = FastaAnnotationParameters()
