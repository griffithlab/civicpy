"""Module for writing CIViC GKS representation to JSON"""

import datetime
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from civicpy.exports.civic_gks_record import (
    CivicGksPredictiveAssertion,
    CivicGksDiagnosticAssertion,
    CivicGksPrognosticAssertion,
)


def get_pkg_version(name: str) -> str:
    """Get package version

    :param name: Name of package to get version for
    :return: Version for package, if found. Otherwise, unknown.
    """
    try:
        return version(name)
    except PackageNotFoundError:
        return "unknown"


class GksOutputMetadata(BaseModel):
    """Define model for GKS JSON Output Metadata"""

    va_spec_python_version: str = Field(
        description="VA-Spec Python version. This can be used to derive the corresponding VA-Spec version.",
        default_factory=lambda: get_pkg_version("ga4gh.va_spec"),
    )
    created_at: str


class GksAssertionError(BaseModel):
    """Define model for representing assertion errors when translating to GKS"""

    assertion_id: int
    message: str


class GksOutput(BaseModel):
    """Define model for representing GKS JSON output"""

    gks_records: list[
        CivicGksPredictiveAssertion
        | CivicGksPrognosticAssertion
        | CivicGksDiagnosticAssertion
    ]
    metadata: GksOutputMetadata
    failed_assertion_ids: list[int] = []
    errors: list[GksAssertionError] = []


class CivicGksWriter:
    """Class for writing CIViC GKS assertions to JSON file

    :param filepath: The output file path to write the JSON file to
    :raises ValueError: If ``filepath`` does not have `.json` suffix
    :param gks_records: List CIViC assertions represented as GKS objects
    """

    def __init__(
        self,
        filepath: Path,
        gks_records: list[
            CivicGksPredictiveAssertion
            | CivicGksPrognosticAssertion
            | CivicGksDiagnosticAssertion
        ],
        errors: list[GksAssertionError] | None = None
    ):
        """Initialize CivicGksWriter class

        :param filepath: The output file path to write the JSON file to
        :raises ValueError: If ``filepath`` does not have `.json` suffix
        :param gks_records: List CIViC assertions represented as GKS objects
        """

        def _default(obj: Any) -> str:
            """Converting datetime objects to ISO-formatted dates (YYYY-MM-DD)

            :param obj: Python object to serialize
            :raises TypeError: If an unsupported type
            :return: JSON string where datetime objects appear as YYYY-MM-DD
            """
            if isinstance(obj, datetime.datetime):
                return obj.isoformat().split("T", 1)[0]

            err_msg = f"Object of type {type(obj)} is not JSON serializable"
            raise TypeError(err_msg)

        if filepath.suffix.lower() != ".json":
            err_msg = "Output file path must end in '.json'."
            raise ValueError(err_msg)

        metadata = GksOutputMetadata(
            created_at=str(datetime.datetime.now(tz=datetime.timezone.utc).strftime("%Y-%m-%d"))
        )

        if errors:
            failed_assertion_ids = [err.assertion_id for err in errors]
        else:
            errors = []
            failed_assertion_ids = []

        output = GksOutput(
            gks_records=gks_records,
            metadata=metadata,
            failed_assertion_ids=failed_assertion_ids,
            errors=errors
        )

        with filepath.open("w+") as wf:
            json.dump(
                output.model_dump(exclude_none=True), wf, indent=2, default=_default
            )
