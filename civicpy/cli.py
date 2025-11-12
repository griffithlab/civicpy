from pathlib import Path
import click
import logging
from civicpy import LOCAL_CACHE_PATH, civic
from civicpy.exports.civic_gks_record import (
    CivicGksRecordError,
    CivicGksPredictiveAssertion,
    CivicGksDiagnosticAssertion,
    CivicGksPrognosticAssertion,
    create_gks_record_from_assertion,
)
from civicpy.exports.civic_gks_writer import CivicGksWriter, GksAssertionError
from civicpy.exports.civic_vcf_writer import CivicVcfWriter
from civicpy.exports.civic_vcf_record import CivicVcfRecord
from civicpy.civic import CoordinateQuery
import vcfpy
from collections import OrderedDict
from civicpy.__version__ import __version__


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
def cli():
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
    "-o",
    "--output-json",
    required=True,
    help="The output file path to write the JSON file to.",
    type=click.Path(
        exists=False,
        readable=True,
        dir_okay=False,
        path_type=Path,
    ),
)
def create_gks_json(organization_id: int, output_json: Path) -> None:
    """Create a JSON file for CIViC assertion records approved by a specific organization that are ready for ClinVar submission, represented as GKS objects.

    For now, we will only support simple molecular profiles and diagnostic, prognostic,
    or predictive assertions.

    :param organization_id: The CIViC organization ID that approved the assertion(s) for submission to ClinVar
    :param output_json: The output file path to write the JSON file to
    """
    try:
        civic.get_organization_by_id(organization_id)
    except Exception:
        logging.exception("Error getting organization %i", organization_id)
        return

    records: list[
        CivicGksDiagnosticAssertion
        | CivicGksPredictiveAssertion
        | CivicGksPrognosticAssertion
    ] = []
    errors: list[GksAssertionError] = []

    for approval in civic.get_all_approvals_ready_for_clinvar_submission_for_org(
        organization_id
    ):
        assertion = approval.assertion
        if assertion.is_valid_for_gks_json(emit_warnings=True):
            try:
                gks_record = create_gks_record_from_assertion(
                    assertion, approval=approval
                )
            except (CivicGksRecordError, NotImplementedError) as e:
                errors.append(
                    GksAssertionError(assertion_id=assertion.id, message=str(e))
                )
                continue
            records.append(gks_record)
        else:
            errors.append(
                GksAssertionError(
                    assertion_id=assertion.id,
                    message="Assertion is not valid for GKS JSON. See logs for more details.",
                )
            )
    if not records:
        logging.warning(
            "No assertions ready for submission to ClinVar found for organization {}".format(
                organization_id
            )
        )
    else:
        CivicGksWriter(output_json, records, errors=errors)


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
