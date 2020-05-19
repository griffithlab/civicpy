import pytest
from civicpy import cli
import tempfile


class TestCli(object):
    def test_create_cache(self):
        tmp_file = tempfile.NamedTemporaryFile('w', delete=False)
        cli.update(['--hard', '--cache-save-path', tmp_file.name], standalone_mode=False)

    def test_create_and_annotate_vcf(self):
        tmp_file = tempfile.NamedTemporaryFile('w', delete=False)
        cli.create_vcf(["-v", tmp_file.name, "--include-status", "accepted"], standalone_mode=False)

        tmp_annotated_file = tempfile.NamedTemporaryFile('w', delete=False)
        cli.annotate_vcf(["--input-vcf", tmp_file.name, "--output-vcf", tmp_annotated_file.name, "--reference", "GRCh37", "--include-status", "accepted"], standalone_mode=False)
