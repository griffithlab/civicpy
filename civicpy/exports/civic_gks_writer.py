"""Module for writing CIViC GKS representation to JSON"""

from datetime import datetime
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


class GksOutput(BaseModel):
    """Define model for representing GKS JSON output"""

    gks_records: list[
        CivicGksPredictiveAssertion
        | CivicGksPrognosticAssertion
        | CivicGksDiagnosticAssertion
    ]
    va_spec_python_version: str = Field(
        description="VA-Spec Python version. This can be used to derive the corresponding VA-Spec version.",
        default_factory=lambda: get_pkg_version("ga4gh.va_spec"),
    )


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

            if isinstance(obj, datetime):
                return obj.isoformat().split("T", 1)[0]

            err_msg = f"Object of type {type(obj)} is not JSON serializable"
            raise TypeError(err_msg)

        if filepath.suffix.lower() != ".json":
            err_msg = "Output file path must end in '.json'."
            raise ValueError(err_msg)

        output = GksOutput(gks_records=gks_records)
        with filepath.open("w+") as wf:
            json.dump(
                output.model_dump(exclude_none=True), wf, indent=2, default=_default
            )
