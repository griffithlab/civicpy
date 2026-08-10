"""Define CIViCpy command-line workflows for cache management and data export."""

from collections import OrderedDict
from collections.abc import Iterable, Iterator
import logging
from pathlib import Path
from typing import TypeAlias

import click
import vcfpy

from civicpy import civic
from civicpy.__env__ import LOCAL_CACHE_PATH
from civicpy.__version__ import __version__
from civicpy.civic import CoordinateQuery
from civicpy.exports.civic_gks_output import GksAssertionError, GksRecord
from civicpy.exports.civic_gks_record import (
    CivicGksRecordError,
    ClinVarSubmissionType,
    create_gks_record_from_assertion,
)
from civicpy.exports.civic_gks_writer import CivicGksWriter
from civicpy.exports.civic_vcf_record import CivicVcfRecord
from civicpy.exports.civic_vcf_writer import CivicVcfWriter
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
_ACCEPTED_STATUS = "accepted"
# An Assertion may have approvals from more than one organization.
_AssertionApprovals: TypeAlias = tuple[civic.Assertion, tuple[civic.Approval, ...]]


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
def cli():
    """Provide CIViCpy command-line tools."""
    pass


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--soft/--hard",
    default=True,
    help="Hard-update from live API (slow) or \
              soft-update from daily precache (fast; default)",
)
@click.option(
    "--cache-save-path",
    help="Filepath to save cache to. Default: {}".format(LOCAL_CACHE_PATH),
    default=LOCAL_CACHE_PATH,
)
def update(soft, cache_save_path):
    """Updates CIViC content from server and stores to local cache file"""
    civic.update_cache(from_remote_cache=soft, local_cache_path=cache_save_path)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "-v", "--vcf-file-path", required=True, help="The file path to write the VCF to."
)
@click.option(
    "-i",
    "--include-status",
    required=True,
    multiple=True,
    type=click.Choice(["accepted", "submitted", "rejected"]),
    help="Limits the variants and annotations in the VCF to only the ones that match the given statuses. \
              May be specified more than once.",
)
def create_vcf(vcf_file_path, include_status):
    """Create a VCF file of CIViC variants"""
    records = []
    for variant in civic.get_all_gene_variants(include_status=include_status):
        if variant.is_valid_for_vcf():
            records.append(CivicVcfRecord(variant, include_status))
    CivicVcfWriter(vcf_file_path, records)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--organization-id",
    required=True,
    help="The CIViC organization ID that approved the assertion(s) for submission to ClinVar.",
    type=int,
)
@click.option(
    "--submission-type",
    type=click.Choice(
        [s.value for s in ClinVarSubmissionType],
        case_sensitive=True,
    ),
    help="The ClinVar submission type to generate GKS JSON for.",
    default=ClinVarSubmissionType.CLINICAL_IMPACT.value,
    show_default=True,
)
@click.option(
    "-o",
    "--output-json",
    required=True,
    help="The output file path to write the dereferenced GKS JSON to.",
    type=click.Path(
        exists=False,
        readable=True,
        dir_okay=False,
        path_type=Path,
    ),
)
def create_gks_json(
    organization_id: int,
    submission_type: str,
    output_json: Path,
) -> None:
    """Export ClinVar-ready Assertions as dereferenced GKS JSON.

    Selects Assertions that the specified CIViC organization approved and marked
    ready for ClinVar submission. Supports simple molecular profiles and diagnostic,
    prognostic, predictive, or oncogenic Assertions. The ClinVar Submission API
    accepts one submission type per request, so ``submission_type`` limits the
    dereferenced output to either clinical impact (diagnostic, prognostic, and
    predictive) or oncogenicity Assertions.

    The Variation Normalizer REST service defaults to
    ``http://127.0.0.1:8000/variation``. Set
    ``CIVICPY_VARIATION_NORMALIZER_URL`` before running the command to use a
    different endpoint.

    ClinVar only supports submitting records of the same submission type for a
    given assertion criteria:
    * Clinical Impact -> diagnostic, prognostic, or predictive assertion
    * Oncogenicity -> oncogenic assertion
    Therefore, create separate GKS JSON files for each submission type.

    :param organization_id: The CIViC organization ID that approved the assertion(s) for submission to ClinVar
    :param submission_type: The ClinVar submission type to generate GKS JSON for.
        Defaults to clinical impact.
    :param output_json: The output file path to write the JSON file to
    """
    if not _organization_exists(organization_id):
        return

    submission_type_filter = ClinVarSubmissionType(submission_type)
    approvals = tuple(
        civic.get_all_approvals_ready_for_clinvar_submission_for_org(organization_id)
    )
    assertions = (approval.assertion for approval in approvals)
    assertion_approvals = _pair_assertions_with_approvals(assertions, approvals)

    _create_gks_export(
        assertion_approvals,
        output_json,
        submission_type=submission_type_filter,
        bundle=False,
    )


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--organization-id",
    required=False,
    help="Only include Assertions approved by this CIViC organization.",
    type=int,
)
@click.option(
    "-o",
    "--output-json",
    required=True,
    help="The output file path to write the referenced GKS bundle to.",
    type=click.Path(
        exists=False,
        readable=True,
        dir_okay=False,
        path_type=Path,
    ),
)
def create_gks_bundle(
    organization_id: int | None,
    output_json: Path,
) -> None:
    """Export all accepted CIViC Assertions as a referenced GKS bundle.

    By default, selects every accepted Assertion and attaches all of its accepted
    Approvals. ``organization_id`` instead limits the bundle to Assertions
    approved by that organization. The bundle contains both clinical significance
    and oncogenicity Statements.

    Set ``CIVICPY_VARIATION_NORMALIZER_URL`` to use a Variation Normalizer
    endpoint other than the default.

    \f
    :param organization_id: Optional CIViC organization whose approved Assertions
        should be included.
    :param output_json: Destination JSON filepath.
    """
    if organization_id is not None and not _organization_exists(organization_id):
        return

    assertions: Iterable[civic.Assertion] = civic.get_all_assertions(
        include_status=[_ACCEPTED_STATUS]
    )
    approvals: Iterable[civic.Approval] = civic.get_all_approvals(
        include_status=[_ACCEPTED_STATUS]
    )

    if organization_id is not None:
        approvals = tuple(
            approval
            for approval in approvals
            if approval.organization_id == organization_id
        )
        approved_assertion_ids = {approval.assertion_id for approval in approvals}
        assertions = (
            assertion
            for assertion in assertions
            if assertion.id in approved_assertion_ids
        )

    assertion_approvals = _pair_assertions_with_approvals(assertions, approvals)

    _create_gks_export(
        assertion_approvals,
        output_json,
        submission_type=None,
        bundle=True,
    )


def _create_gks_export(
    assertion_approvals: Iterable[_AssertionApprovals],
    output_json: Path,
    submission_type: ClinVarSubmissionType | None,
    bundle: bool,
) -> None:
    """Transform eligible CIViC Assertions into GKS Statements and write JSON.

    The Statements are written either as the default dereferenced document or as
    a referenced CIViC GKS bundle, according to ``bundle``.

    :param assertion_approvals: Assertions paired with their applicable approvals.
    :param output_json: Destination JSON filepath.
    :param submission_type: ClinVar submission filter, or ``None`` to include all
        supported Statement types.
    :param bundle: If ``True``, write a referenced CIViC GKS bundle; otherwise
        write the default dereferenced GKS JSON document.
    """
    gks_records, errors = _transform_assertions_to_gks(
        assertion_approvals,
        submission_type,
    )

    if not gks_records:
        logging.warning("No eligible Assertions found for GKS export")
        return

    CivicGksWriter(output_json, gks_records, errors=errors, bundle=bundle)


def _transform_assertions_to_gks(
    assertion_approvals: Iterable[_AssertionApprovals],
    submission_type: ClinVarSubmissionType | None,
) -> tuple[list[GksRecord], list[GksAssertionError]]:
    """Transform eligible CIViC Assertions into VA-Spec GKS Statement models.

    :param assertion_approvals: Assertions paired with their approvals.
    :param submission_type: Optional ClinVar submission-type filter.
    :return: Successfully transformed GKS Statements and errors for Assertions
        that could not be transformed.
    """
    gks_records: list[GksRecord] = []
    errors: list[GksAssertionError] = []

    for assertion, approvals in assertion_approvals:
        if not assertion.is_valid_for_gks_json(emit_warnings=True):
            errors.append(
                GksAssertionError(
                    assertion_id=assertion.id,
                    message="Assertion is not valid for GKS JSON. See logs for more details.",
                )
            )
            continue

        try:
            gks_record = create_gks_record_from_assertion(
                assertion,
                approval=approvals,
                submission_type_filter=submission_type,
            )
        except (CivicGksRecordError, NotImplementedError) as error:
            errors.append(
                GksAssertionError(assertion_id=assertion.id, message=str(error))
            )
            continue

        gks_records.append(gks_record)

    return gks_records, errors


def _pair_assertions_with_approvals(
    assertions: Iterable[civic.Assertion],
    approvals: Iterable[civic.Approval],
) -> Iterator[_AssertionApprovals]:
    """Pair each unique Assertion with all of its approvals.

    Assertions without approvals are retained. Assertions and approvals are
    sorted to make output deterministic regardless of API response ordering.

    :param assertions: Assertions to include.
    :param approvals: CIViC Approvals to group.
    :return: Assertions paired with their sorted approvals.
    """
    approvals_by_assertion: dict[int, list[civic.Approval]] = {}
    for approval in approvals:
        assertion_id = approval.assertion_id
        approvals_by_assertion.setdefault(assertion_id, []).append(approval)

    assertions_by_id = {assertion.id: assertion for assertion in assertions}
    for assertion_id in sorted(assertions_by_id):
        assertion_approvals = tuple(
            sorted(
                approvals_by_assertion.get(assertion_id, []),
                key=lambda approval: (approval.organization_id, approval.id),
            )
        )
        yield assertions_by_id[assertion_id], assertion_approvals


def _organization_exists(organization_id: int) -> bool:
    """Return whether a CIViC organization can be retrieved.

    :param organization_id: CIViC organization identifier to validate.
    :return: ``True`` when the organization exists; otherwise ``False``.
    """
    try:
        civic.get_organization_by_id(organization_id)
    except Exception:
        logging.exception("Error getting organization %i", organization_id)
        return False
    return True


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--input-vcf", required=True, help="A VCF to annotate with information from CIViC."
)
@click.option(
    "--output-vcf", required=True, help="The file path to write the annotated VCF to."
)
@click.option(
    "--reference",
    required=True,
    type=click.Choice(["NCBI36", "GRCh37", "GRCh38"]),
    help="The reference sequence build used to create the input VCF",
)
@click.option(
    "-i",
    "--include-status",
    required=True,
    multiple=True,
    type=click.Choice(["accepted", "submitted", "rejected"]),
    help="Limits the variants and annotations in the VCF to only the ones that match the given statuses. \
              May be specified more than once.",
)
def annotate_vcf(input_vcf, output_vcf, reference, include_status):
    """Annotate a VCF with information from CIViC"""
    reader = vcfpy.Reader.from_path(input_vcf)
    new_header = reader.header.copy()
    new_header.add_info_line(
        OrderedDict(
            [
                ("ID", "CIVIC"),
                ("Number", "."),
                ("Type", "String"),
                ("Description", CivicVcfWriter.CSQ_DESCRIPTION),
            ]
        )
    )
    writer = vcfpy.Writer.from_path(output_vcf, new_header)
    for entry in reader:
        for alt in entry.ALT:
            position = entry.POS
            ref = entry.REF
            alt = alt.value
            if len(ref) == 1 and len(alt) == 1:
                start = position
                end = position
            else:
                if len(ref) == len(alt):
                    start = position
                    end = position + len(ref) - 1
                else:
                    alt = alt[1:]
                    ref = ref[1:]
                    if len(ref) > len(alt):
                        start = position + 1
                        end = start + len(ref) - 1
                        if alt == "":
                            alt = None
                    else:
                        start = position
                        if ref == "":
                            ref = None
                            end = start + 1
                        else:
                            end = start + len(ref) - 1
            query = CoordinateQuery(entry.CHROM, start, end, alt, ref, reference)
            variants = civic.search_variants_by_coordinates(query, search_mode="exact")
            if variants is not None:
                if len(variants) == 1:
                    record = CivicVcfRecord(variants[0], include_status)
                    csq = record.INFO["CSQ"]
                    if len(csq) > 0:
                        entry.INFO["CIVIC"] = csq
                elif len(variants) > 1:
                    print(
                        "More than one variant found for start {} stop {} ref {} alt {}. CIViC Variants IDs: {}".format(
                            start,
                            end,
                            ref,
                            alt,
                            ",".join(list(map(lambda v: str(v.id), variants))),
                        )
                    )
            writer.write_record(entry)
    writer.close()
    reader.close()


if __name__ == "__main__":
    cli()
