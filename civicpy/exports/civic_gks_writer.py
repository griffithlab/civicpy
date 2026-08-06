"""Write CIViC GKS assertions as dereferenced or referenced JSON.

This module owns the dereferenced output model and the file-writing entry point.
It delegates referenced document construction to
:mod:`civicpy.exports.civic_gks_bundle_builder` and uses the bundle models from
:mod:`civicpy.exports.civic_gks_bundle`.
"""

import datetime
import json
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from civicpy.exports.civic_gks_bundle import (
    GksBundleOutput,
)
from civicpy.exports.civic_gks_bundle_builder import build_gks_bundle
from civicpy.exports.civic_gks_output import (
    GksAssertionError,
    GksOutputMetadata,
    GksRecord,
)


def _serialize_json_default(value: Any) -> str:
    """Serialize datetime values consistently in GKS JSON documents.

    :param value: Value passed to :func:`json.dump` that needs serialization.
    :raises TypeError: If ``value`` is not a supported datetime object.
    :return: An ISO-8601 calendar date.
    """
    if isinstance(value, datetime.datetime):
        return value.date().isoformat()
    raise TypeError(f"Object of type {type(value)} is not JSON serializable")


class GksOutput(BaseModel):
    """Dereferenced export with each assertion's related objects inline."""

    gks_records: list[GksRecord]
    metadata: GksOutputMetadata
    failed_assertion_ids: list[int]
    errors: list[GksAssertionError]


class CivicGksWriter:
    """Write validated CIViC GKS assertions to a JSON document.

    ``bundle=False`` writes the default dereferenced format, with related objects
    inlined in each assertion. ``bundle=True`` references shared GKS objects from
    keyed collections with JSON Pointers.
    """

    def __init__(
        self,
        filepath: Path,
        gks_records: list[GksRecord],
        errors: list[GksAssertionError] | None = None,
        bundle: bool = False,
    ) -> None:
        """Write CIViC GKS Statements to a dereferenced or bundled JSON file.

        :param filepath: Destination JSON filepath.
        :param gks_records: VA-Spec GKS Statements translated from CIViC
            Assertions.
        :param errors: Assertions that could not be represented in the export.
        :param bundle: If ``True``, write the referenced GKS Bundle Format;
            otherwise write the default dereferenced format.
        :raises ValueError: If ``filepath`` does not use the ``.json`` suffix.
        """
        if filepath.suffix.lower() != ".json":
            raise ValueError("Output file path must end in '.json'.")

        metadata = GksOutputMetadata(
            created_at=datetime.datetime.now(tz=datetime.timezone.utc).strftime(
                "%Y-%m-%d"
            )
        )
        export_errors = errors or []

        if bundle:
            output: GksBundleOutput | GksOutput = build_gks_bundle(
                gks_records, metadata, export_errors
            )
        else:
            output = GksOutput(
                gks_records=gks_records,
                metadata=metadata,
                failed_assertion_ids=[error.assertion_id for error in export_errors],
                errors=export_errors,
            )

        with filepath.open("w", encoding="utf-8") as write_file:
            json.dump(
                output.model_dump(exclude_none=True),
                write_file,
                indent=2,
                default=_serialize_json_default,
            )
