"""Module to support normalizing molecular profiles to VRS objects via the VICC
Variation Normalizer (https://github.com/cancervariants/variation-normalization/).

Only a REST data proxy is provided by default to avoid adding additional dependencies and
external resource setup that the Variation Normalizer Python API requires. Applications
may implement their own Variation Normalizer Python data proxy.
"""

import logging
import re
from abc import ABC, abstractmethod
from os import getenv
from urllib.parse import unquote

import requests
from ga4gh.vrs.models import Allele, CopyNumberChange, VrsType
from pydantic import ValidationError
from pydantic.dataclasses import dataclass

from civicpy.civic import GeneVariant, MolecularProfile

_logger = logging.getLogger(__name__)

MP_NAME_PATTERN = (
    r"(?P<gene>[\w-]+)"
    r"(?:\s+(?P<change>[^\s]+))?"
    r"(?:\s+(?P<c_change>\(?c\..*\)?))?"
)

# Variant names that are known to not be supported in the Variation Normalizer
VAR_NAMES_NOT_SUPPORTED = {
    "mutation",
    "exon",
    "overexpression",
    "frameshift",
    "promoter",
    "deletion",
    "type",
    "insertion",
    "expression",
    "duplication",
    "copy",
    "underexpression",
    "number",
    "variation",
    "repeat",
    "rearrangement",
    "activation",
    "mislocalization",
    "translocation",
    "wild",
    "polymorphism",
    "frame",
    "shift",
    "loss",
    "function",
    "levels",
    "inactivation",
    "snp",
    "fusion",
    "dup",
    "truncation",
    "homozygosity",
    "gain",
    "phosphorylation",
}


@dataclass
class MolecularProfileNameComponents:
    """Parsed molecular profile name components."""

    gene: str
    change: str | None
    c_change: str | None


class VariationNormalizerError(Exception):
    """Raised when variation normalization fails."""


class VariationNormalizerDataProxy(ABC):
    """Abstract base class for the Variation Normalizer."""

    @abstractmethod
    def normalize(
        self,
        expr: str,
    ) -> Allele | CopyNumberChange | None:
        """Normalize a variation expression via the Variation Normalizer

        :param expr: Variation expression to normalize.
        :return: Normalized VRS variation object if normalization succeeds.
            Otherwise ``None``.
        :raises VariationNormalizerError: If the request fails or the response is not a
            valid VRS variation object.
        """
        raise NotImplementedError

    @staticmethod
    def _extract_mp_components(
        molecular_profile_name: str,
    ) -> MolecularProfileNameComponents | None:
        """Extract components from a molecular profile name.

        :param molecular_profile_name: CIViC Molecular Profile name.
        :return: Molecular profile name components if the pattern matches.  Otherwise
            ``None``.
        """
        match = re.fullmatch(
            MP_NAME_PATTERN,
            molecular_profile_name.strip(),
        )

        if not match:
            return None

        return MolecularProfileNameComponents(**match.groupdict())

    @staticmethod
    def _is_supported_change(
        change: str,
        molecular_profile_id: int,
    ) -> bool:
        """Determine whether a change is supported by the Variation Normalizer.

        :param change: Variant change component extracted from the molecular
            profile name. This is what will be passed to the Variation Normalizer.
        :param molecular_profile_id: CIViC Molecular Profile ID. Used for logging.
        :return: ``True`` if the change appears to be supported by the
            Variation Normalizer. Otherwise ``False``.
        """
        change_lower = change.lower()

        if (
            change_lower.endswith("fs")
            or any(c in change_lower for c in ("-", "/"))
            or bool(set(change_lower.split()) & VAR_NAMES_NOT_SUPPORTED)
        ):
            _logger.warning(
                "Unsupported molecular profile change for Variation Normalizer. mpid=%i, change='%s'",
                molecular_profile_id,
                change,
            )
            return False

        return True

    def _build_normalizable_query(
        self,
        molecular_profile: MolecularProfile,
    ) -> str | None:
        """Build a Variation Normalizer query (of the form ``{gene} {change}``) from a
        molecular profile.

        :param molecular_profile: CIViC Molecular Profile.
        :return: Query use in the Variation Normalizer if supported. Otherwise
            ``None``
        """
        variants = molecular_profile.variants
        mp_id = molecular_profile.id
        mp_name = molecular_profile.name

        if not variants:
            _logger.warning("No variants found. mpid=%i, name='%s'", mp_id, mp_name)
            return None

        if len(variants) > 1:
            _logger.warning(
                "Complex molecular profiles are not supported. mpid=%i, name='%s'",
                mp_id,
                mp_name,
            )
            return None

        if not isinstance(
            variants[0],
            GeneVariant,
        ):
            _logger.warning(
                "Variant type '%s' is not supported. mpid=%i, name='%s'",
                variants[0].__class__.__name__,
                mp_id,
                mp_name,
            )
            return None

        components = self._extract_mp_components(mp_name)

        if not components:
            _logger.warning(
                "Unable to parse molecular profile name. mpid=%i, name='%s'",
                mp_id,
                mp_name,
            )
            return None

        if components.c_change:
            _logger.warning(
                "Molecular profiles containing cDNA changes are not supported. mpid=%i, name='%s'",
                mp_id,
                mp_name,
            )
            return None

        change = components.change

        if not change:
            _logger.warning(
                "No change component found in molecular profile name. mpid=%i, name='%s'",
                mp_id,
                mp_name,
            )
            return None

        if not self._is_supported_change(
            change,
            molecular_profile.id,
        ):
            return None

        return f"{components.gene} {change}"

    def normalize_molecular_profile(
        self,
        molecular_profile: MolecularProfile,
    ) -> Allele | CopyNumberChange | None:
        """Normalize a CIViC Molecular Profile to a VRS variation.

        :param molecular_profile: CIViC Molecular Profile.
        :return: Normalized VRS variation if successful. Otherwise ``None``.
        """
        expression = self._build_normalizable_query(molecular_profile)
        if not expression:
            return None

        normalized_variation = self.normalize(expression)

        if not normalized_variation:
            _logger.warning(
                "Variation Normalizer failed to normalize query. mpid=%i, query='%s'",
                molecular_profile.id,
                expression,
            )
            return None

        return normalized_variation


class VariationNormalizerRestDataProxy(VariationNormalizerDataProxy):
    """REST data proxy for the Variation Normalizer service."""

    DEFAULT_BASE_URL = "http://127.0.0.1:8000/variation"
    BASE_URL_ENV_VAR = "CIVICPY_VARIATION_NORMALIZER_URL"

    def __init__(self, base_url: str | None = None):
        """Initialize the Variation Normalizer REST client.

        :param base_url: Base URL for the Variation Normalizer REST service.
            If not provided, ``CIVICPY_VARIATION_NORMALIZER_URL`` environment
            variable will be used, followed by the ``DEFAULT_BASE_URL``.
        """
        self.base_url = (
            base_url or getenv(self.BASE_URL_ENV_VAR) or self.DEFAULT_BASE_URL
        )

    @staticmethod
    def _validate_vrs_variation(
        variation: dict,
        query: str,
    ) -> Allele | CopyNumberChange | None:
        """Validate and construct a VRS variation object.

        :param variation: Variation object returned by the Variation Normalizer.
        :param query: Original normalization query. Used for logging.
        :return: Validated VRS variation object if supported. Otherwise ``None``.
        :raise VariationNormalizerError: If the variation cannot be validated as a
            supported VRS variation object.
        """
        variation_type = variation.get("type")

        try:
            if variation_type == VrsType.ALLELE:
                return Allele.model_validate(variation)

            if variation_type == VrsType.CN_CHANGE:
                return CopyNumberChange.model_validate(variation)

            return None

        except ValidationError as e:
            msg = f"Variation Normalizer returned invalid VRS variation object. query={query!r}, variation_type={variation_type!r}"
            _logger.exception(msg)
            raise VariationNormalizerError(msg) from e

    def normalize(
        self,
        expr: str,
    ) -> Allele | CopyNumberChange | None:
        """Normalize a variation expression via the ``/normalize`` endpoint

        :param expr: Variation expression to normalize.
        :return: Normalized VRS variation object if normalization succeeds.
            Otherwise ``None``.
        :raises VariationNormalizerError: If the request fails or the response is not a
            valid VRS variation object.
        """
        try:
            response = requests.get(
                f"{self.base_url}/normalize",
                params={"q": unquote(expr)},
                timeout=15,
            )
            response.raise_for_status()

        except requests.exceptions.HTTPError as e:
            msg = (
                "Variation Normalizer returned unexpected HTTP status. "
                f"query={expr!r}, "
                f"status_code={e.response.status_code}"
            )
            _logger.exception(msg)
            raise VariationNormalizerError(msg) from e

        except requests.exceptions.RequestException as e:
            msg = f"Variation Normalizer request failed. query={expr!r}, error={e!s}"
            _logger.exception(msg)
            raise VariationNormalizerError(msg) from e

        data = response.json()

        variation = data.get("variation")
        if not variation:
            _logger.warning(
                "Variation Normalizer returned no variation object. query=%r",
                expr,
            )
            return None

        variation["name"] = expr
        variation.pop("extensions", None)

        return self._validate_vrs_variation(
            variation,
            expr,
        )
