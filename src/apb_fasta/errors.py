"""Expected APB FASTA annotation failures."""


class FastaAnnotationError(ValueError):
    """The supplied APB2 result or protein frame cannot be annotated safely."""
