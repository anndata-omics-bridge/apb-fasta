"""User-selected FASTA annotation behavior."""

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True, slots=True)
class FastaAnnotationParameters:
    """Configuration shared by peptide validation and protein annotation."""

    matcher_backend: Literal["auto", "ahocorapy", "ahocorasick_rs"] = "auto"


DEFAULT_FASTA_ANNOTATION_PARAMETERS = FastaAnnotationParameters()
