"""Define data shared by dereferenced and bundled CIViC GKS outputs.

This module owns the common record type, metadata, and error models. It does not
construct bundle references or write files; those responsibilities belong to
:mod:`civicpy.exports.civic_gks_bundle` and
:mod:`civicpy.exports.civic_gks_writer`, respectively.
"""

from importlib.metadata import PackageNotFoundError, version
from typing import TypeAlias

from pydantic import BaseModel, Field

from civicpy.exports.civic_gks_record import (
    CivicGksClinSigAssertion,
    CivicGksOncogenicAssertion,
)

GksRecord: TypeAlias = CivicGksClinSigAssertion | CivicGksOncogenicAssertion


def _get_pkg_version(name: str) -> str:
    """Return an installed package version, or ``unknown`` when unavailable.

    :param name: Distribution name to look up.
    :return: Installed distribution version, or ``unknown``.
    """
    try:
        return version(name)
    except PackageNotFoundError:
        return "unknown"


class GksOutputMetadata(BaseModel):
    """Generation date and VA-Spec implementation version for a GKS export."""

    va_spec_python_version: str = Field(
        description="VA-Spec Python version used to generate this export.",
        default_factory=lambda: _get_pkg_version("ga4gh.va_spec"),
    )
    created_at: str


class GksAssertionError(BaseModel):
    """Describe a CIViC Assertion that could not become a GKS Statement."""

    assertion_id: int
    message: str
