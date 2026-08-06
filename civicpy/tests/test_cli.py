from datetime import datetime
from pathlib import Path
from typing import Any
import pytest
from unittest.mock import Mock, patch
from civicpy import cli, civic
import tempfile
import json

from civicpy.exports.civic_gks_bundle import GksBundleOutput
from civicpy.exports.civic_gks_writer import GksOutput


def _bundle_assertion_ids(bundle: dict[str, Any]) -> set[str]:
    """Return the CIViC Assertion IDs stored in a bundle."""
    return {
        identifier
        for identifier in bundle["statement"]
        if identifier.startswith("civic.aid:")
    }


def check_metadata(metadata: dict[str, Any], bundle: bool = False) -> None:
    """Check that metadata output is correct"""
    expected_keys = {"va_spec_python_version", "created_at"}
    if bundle:
        expected_keys.update(
            {
                "bundle_format",
                "bundle_format_version",
                "statistics",
            }
        )
        assert metadata["bundle_format"] == "civic-gks-bundle"
        assert metadata["bundle_format_version"] == "0.1.0"
        statistics = metadata["statistics"]
        assert set(statistics) == {"collections"}
        assert all(
            collection["count"] >= 0
            and all(count >= 0 for count in collection.get("types", {}).values())
            for collection in statistics["collections"].values()
        )
    assert set(metadata.keys()) == expected_keys
    va_spec_python_version = metadata["va_spec_python_version"]
    assert isinstance(va_spec_python_version, str) and va_spec_python_version

    created_at = metadata["created_at"]
    assert datetime.strptime(created_at, "%Y-%m-%d")


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

    @patch("civicpy.civic.get_all_assertions")
    @patch("civicpy.civic.get_all_approvals")
    @patch("civicpy.civic.get_all_approvals_ready_for_clinvar_submission_for_org")
    @patch("civicpy.civic.get_assertion_by_id", wraps=civic.get_assertion_by_id)
    @pytest.mark.parametrize("bundle", [False, True])
    def test_create_gks_exports_assertions_found(
        self,
        mock_assertion,
        mock_clinvar_approvals,
        mock_all_approvals,
        mock_all_assertions,
        mocked_normalizer,
        bundle,
    ):
        """Test dereferenced and bundle CLI exports for eligible Assertions."""
        assertions = {
            assertion_id: civic.get_assertion_by_id(assertion_id)
            for assertion_id in (6, 202)
        }
        mock_assertion.side_effect = assertions.__getitem__
        mock_all_assertions.return_value = list(assertions.values())
        approvals = [
            civic.Approval(
                type="approval",
                id=approval_id,
                status="ACTIVE",
                ready_for_clinvar_submission=True,
                organization_id=1,
                last_reviewed="2025-03-19T13:38:54Z",
                assertion_id=assertion_id,
                clinvar_accession=("SCV000000001" if assertion_id == 6 else None),
                partial=True,
            )
            for approval_id, assertion_id in ((1, 6), (2, 202))
        ]
        mock_all_approvals.return_value = [
            approvals[0],
            civic.Approval(
                type="approval",
                id=3,
                status="ACTIVE",
                ready_for_clinvar_submission=False,
                organization_id=2,
                last_reviewed="2025-04-20T10:30:00Z",
                assertion_id=6,
                clinvar_accession="SCV000000002",
                partial=True,
            ),
        ]
        mock_clinvar_approvals.return_value = approvals

        with tempfile.NamedTemporaryFile("w+", suffix=".json", delete=True) as tmp_file:
            try:
                with patch(
                    "civicpy.cli.VariationNormalizerRESTDataProxy",
                    return_value=mocked_normalizer,
                ):
                    command = ["-o", Path(tmp_file.name)]
                    if not bundle:
                        command[0:0] = [
                            "--organization-id",
                            1,
                            "--submission-type",
                            "oncogenicity",
                        ]
                    export_command = (
                        cli.create_gks_bundle if bundle else cli.create_gks_json
                    )
                    export_command(command)
            except SystemExit as e:
                assert e.code == 0

            with open(tmp_file.name, "r") as f:
                gks_output = json.load(f)
                expected_model = GksBundleOutput if bundle else GksOutput
                assert set(gks_output.keys()) == set(expected_model.model_fields.keys())
                check_metadata(gks_output["metadata"], bundle=bundle)
                if bundle:
                    GksBundleOutput.model_validate(gks_output)
                    statement = gks_output["statement"]["civic.aid:6"]
                    accessions_by_contributor = {
                        contribution["contributor"]: contribution["extensions"][0]
                        for contribution in statement["contributions"]
                    }
                    assert _bundle_assertion_ids(gks_output) == {
                        "civic.aid:6",
                        "civic.aid:202",
                    }
                    assert set(gks_output["agent"]) == {
                        "civic.organization:1",
                        "civic.organization:2",
                    }
                    assert statement["extensions"] == [
                        {
                            "name": "clinvarAccessions",
                            "value": ["SCV000000001", "SCV000000002"],
                        }
                    ]
                    assert accessions_by_contributor == {
                        "#/agent/civic.organization:1": {
                            "name": "clinvarAccession",
                            "value": "SCV000000001",
                        },
                        "#/agent/civic.organization:2": {
                            "name": "clinvarAccession",
                            "value": "SCV000000002",
                        },
                    }
                else:
                    assert gks_output["failed_assertion_ids"] == [6]
                    assert [
                        record["id"] for record in gks_output["gks_records"]
                    ] == ["civic.aid:202"]

    @patch("civicpy.civic.get_all_assertions")
    @patch("civicpy.civic.get_all_approvals")
    @patch("civicpy.civic.get_assertion_by_id", wraps=civic.get_assertion_by_id)
    def test_create_gks_bundle_filters_by_organization(
        self,
        mock_assertion: Mock,
        mock_approvals: Mock,
        mock_all_assertions: Mock,
        mocked_normalizer: Mock,
        tmp_path: Path,
    ) -> None:
        """Test the optional organization filter for bundle output."""
        assertions = {
            assertion_id: civic.get_assertion_by_id(assertion_id)
            for assertion_id in (6, 202)
        }
        mock_assertion.side_effect = assertions.__getitem__
        mock_all_assertions.return_value = list(assertions.values())
        mock_approvals.return_value = [
            civic.Approval(
                type="approval",
                id=approval_id,
                status="ACTIVE",
                ready_for_clinvar_submission=False,
                organization_id=organization_id,
                last_reviewed="2025-03-19T13:38:54Z",
                assertion_id=assertion_id,
                partial=True,
            )
            for approval_id, organization_id, assertion_id in (
                (1, 1, 6),
                (2, 2, 202),
            )
        ]
        output_path = tmp_path / "civic-gks-bundle.json"

        with patch(
            "civicpy.cli.VariationNormalizerRESTDataProxy",
            return_value=mocked_normalizer,
        ):
            try:
                cli.create_gks_bundle(["--organization-id", 1, "-o", output_path])
            except SystemExit as error:
                assert error.code == 0

        with output_path.open() as read_file:
            bundle = json.load(read_file)
        assert _bundle_assertion_ids(bundle) == {"civic.aid:6"}
        assert set(bundle["agent"]) == {"civic.organization:1"}

    @patch("civicpy.civic.get_all_approvals_ready_for_clinvar_submission_for_org")
    def test_create_gks_json_assertions_not_valid(
        self, mock_approvals, mocked_normalizer
    ):
        """Test that CLI create_gks_json works as expected when assertion is not valid for GKS JSON"""
        mock_approvals.return_value = [
            civic.Approval(
                type="approval",
                id=1,
                status="ACTIVE",
                ready_for_clinvar_submission=True,
                organization_id=1,
                last_reviewed="2025-03-19T13:38:54Z",
                assertion_id=4,  # Assertion is predisposing, which is not supported
                partial=True,
            ),
            civic.Approval(
                type="approval",
                id=2,
                status="ACTIVE",
                ready_for_clinvar_submission=True,
                organization_id=1,
                last_reviewed="2025-03-19T13:38:54Z",
                assertion_id=6,
                partial=True,
            ),
        ]

        with tempfile.NamedTemporaryFile("w+", suffix=".json", delete=True) as tmp_file:
            try:
                with patch(
                    "civicpy.cli.VariationNormalizerRESTDataProxy",
                    return_value=mocked_normalizer,
                ):
                    cli.create_gks_json(
                        ["--organization-id", 1, "-o", Path(tmp_file.name)]
                    )
            except SystemExit as e:
                assert e.code == 0

            with open(tmp_file.name, "r") as f:
                gks_output = json.load(f)
                check_metadata(gks_output["metadata"])
                assert [
                    record["id"] for record in gks_output["gks_records"]
                ] == ["civic.aid:6"]
                assert gks_output["failed_assertion_ids"] == [4]
                assert gks_output["errors"] == [
                    {
                        "assertion_id": 4,
                        "message": "Assertion is not valid for GKS JSON. See logs for more details.",
                    }
                ]

    def test_create_gks_json_no_organization(self, tmp_path, caplog):
        """Test that CLI create_gks_json works as expected when organization ID does not exist"""
        output_file = tmp_path / "gks.json"

        try:
            cli.create_gks_json(["--organization-id", 99999999, "-o", output_file])
        except SystemExit as e:
            assert e.code == 0

        assert not output_file.exists()
        assert "Error getting organization 99999999" in caplog.text

    @patch("civicpy.cli.VariationNormalizerRESTDataProxy")
    @patch("civicpy.civic.get_all_approvals_ready_for_clinvar_submission_for_org")
    def test_create_gks_json_variation_normalizer_url(
        self, mock_approvals, mock_normalizer, tmp_path
    ):
        """The CLI configures one normalizer proxy for the export operation."""
        mock_approvals.return_value = []
        output_file = tmp_path / "gks.json"
        normalizer_url = "http://variation-normalizer.example/variation"

        try:
            cli.create_gks_json(
                [
                    "--organization-id",
                    1,
                    "--variation-normalizer-url",
                    normalizer_url,
                    "-o",
                    output_file,
                ]
            )
        except SystemExit as e:
            assert e.code == 0

        mock_normalizer.assert_called_once_with(normalizer_url)

    @patch("civicpy.civic.get_all_approvals_ready_for_clinvar_submission_for_org")
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
        assert "No eligible Assertions found for GKS export" in caplog.text
