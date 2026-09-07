"""Package surface tests."""


def test_package_imports() -> None:
    import apb_fasta.annotation
    import apb_fasta.cli
    import apb_fasta.configuration

    assert apb_fasta.annotation.FastaAnnotationParser
    assert apb_fasta.cli.app
    assert apb_fasta.configuration.FastaAnnotationParameters
