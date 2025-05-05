from pathlib import Path
import pytest
from unittest.mock import patch
from civicpy import cli, civic
import tempfile
import json


class TestCli(object):
    @pytest.mark.skip(reason="Long running test")
    def test_create_cache(self):
        tmp_file = tempfile.NamedTemporaryFile("w", delete=False)
        cli.update(
            ["--hard", "--cache-save-path", tmp_file.name], standalone_mode=False
        )

    @pytest.mark.skip(reason="Long running test")
    def test_create_and_annotate_vcf(self):
        tmp_file = tempfile.NamedTemporaryFile("w", delete=False)
        cli.create_vcf(
            ["-v", tmp_file.name, "--include-status", "accepted"], standalone_mode=False
        )

    @patch("civicpy.civic.get_all_endorsements_ready_for_clinvar_submission_for_org")
    @patch("civicpy.civic.get_assertion_by_id", wraps=civic.get_assertion_by_id)
    def test_create_gks_json_assertions_found(self, mock_assertion, mock_endorsements):
        """Test that CLI create_gks_json works as expected when assertions are ready for clinvar submission"""
        mock_assertion.return_value = civic.get_assertion_by_id(6)
        mock_endorsements.return_value = [
            civic.Endorsement(
                type="endorsement",
                id=1,
                status="ACTIVE",
                ready_for_clinvar_submission=True,
                organization_id=1,
                last_reviewed="2025-03-19T13:38:54Z",
                assertion_id=6,
                partial=True,
            )
        ]

        with tempfile.NamedTemporaryFile("w+", suffix=".json", delete=True) as tmp_file:
            try:
                cli.create_gks_json(["--organization-id", 1, "-o", Path(tmp_file.name)])
            except SystemExit as e:
                assert e.code == 0

            with open(tmp_file.name, "r") as f:
                gks_output = json.load(f)
                assert set(gks_output.keys()) == {
                    "gks_records",
                    "va_spec_python_version",
                }

                assert len(gks_output["gks_records"]) == 1
                assert gks_output["gks_records"][0]["id"] == "civic.aid:6"
                va_spec_python_version = gks_output["va_spec_python_version"]
                assert (
                    isinstance(va_spec_python_version, str) and va_spec_python_version
                )

    def test_create_gks_json_no_organization(
        self, tmp_path, caplog
    ):
        """Test that CLI create_gks_json works as expected when organization ID does not exist"""
        output_file = tmp_path / "gks.json"

        try:
            cli.create_gks_json(["--organization-id", 99999999, "-o", output_file])
        except SystemExit as e:
            assert e.code == 0

        assert not output_file.exists()
        assert (
            "Error getting organization 99999999"
            in caplog.text
        )

    @patch("civicpy.civic.get_all_endorsements_ready_for_clinvar_submission_for_org")
    def test_create_gks_json_no_assertions_found(
        self, mock_assertions, tmp_path, caplog
    ):
        """Test that CLI create_gks_json works as expected when assertions are not ready for clinvar submission"""
        mock_assertions.return_value = []
        output_file = tmp_path / "gks.json"

        try:
            cli.create_gks_json(["--organization-id", 1, "-o", output_file])
        except SystemExit as e:
            assert e.code == 0

        assert not output_file.exists()
        assert (
            "No assertions ready for submission to ClinVar found for organization 1"
            in caplog.text
        )
