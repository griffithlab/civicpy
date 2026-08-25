"""Define data shared by dereferenced and bundled CIViC GKS outputs.

This module owns the common record type, metadata, and error models. It does not
construct bundle references or write files; those responsibilities belong to
:mod:`civicpy.exports.gks.bundle` and
:mod:`civicpy.exports.civic_gks_writer`, respectively.
"""

from importlib.metadata import PackageNotFoundError, version
from typing import TypeAlias

from ga4gh.cat_vrs import CATVRS_VERSION
from ga4gh.core import CORE_VERSION
from ga4gh.va_spec import VASPEC_VERSION
from ga4gh.vrs import VRS_VERSION
from pydantic import BaseModel, ConfigDict, Field
from pydantic.alias_generators import to_camel

from civicpy.exports.civic_gks_record import (
    CivicGksClinSigAssertion,
    CivicGksOncogenicAssertion,
)

GksRecord: TypeAlias = CivicGksClinSigAssertion | CivicGksOncogenicAssertion


class GksModel(BaseModel):
    """Base model using camelCase aliases for GKS JSON fields."""

    model_config = ConfigDict(alias_generator=to_camel, populate_by_name=True)


def _get_pkg_version(name: str) -> str:
    """Return an installed package version, or ``unknown`` when unavailable.

    :param name: Distribution name to look up.
    :return: Installed distribution version, or ``unknown``.
    """
    try:
        return version(name)
    except PackageNotFoundError:
        return "unknown"


class ImplementationVersions(GksModel):
    """Python package versions used to generate a GKS export."""

    vrs_python: str = Field(
        alias="VRSPython",
        description="VRS-Python version used to generate this export.",
        default_factory=lambda: _get_pkg_version("ga4gh.vrs"),
    )
    cat_vrs_python: str = Field(
        alias="CatVRSPython",
        description="Cat-VRS Python version used to generate this export.",
        default_factory=lambda: _get_pkg_version("ga4gh.cat_vrs"),
    )
    va_spec_python: str = Field(
        alias="VASpecPython",
        description="VA-Spec Python version used to generate this export.",
        default_factory=lambda: _get_pkg_version("ga4gh.va_spec"),
    )


class SpecificationVersions(GksModel):
    """GA4GH specification versions represented in a GKS export."""

    gks_core: str = Field(
        alias="GKSCore",
        description="GKS-Core version represented in this export.",
        default=CORE_VERSION,
    )
    vrs: str = Field(
        alias="VRS",
        description="VRS version represented in this export.",
        default=VRS_VERSION,
    )
    cat_vrs: str = Field(
        alias="CatVRS",
        description="Cat-VRS version represented in this export.",
        default=CATVRS_VERSION,
    )
    va_spec: str = Field(
        alias="VASpec",
        description="VA-Spec version used to generate this export.",
        default=VASPEC_VERSION,
    )


class GksOutputMetadata(GksModel):
    """Generation date and version provenance for a GKS export."""

    implementation_versions: ImplementationVersions = Field(
        alias="implementationVersions",
        default_factory=ImplementationVersions,
    )
    specification_versions: SpecificationVersions = Field(
        alias="specificationVersions",
        default_factory=SpecificationVersions,
    )
    created_at: str


class GksAssertionError(GksModel):
    """Describe a CIViC Assertion that could not become a GKS Statement."""

    assertion_id: int
    message: str
