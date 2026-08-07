"""Normalize molecular profiles to VRS objects with the VICC Variation Normalizer.

The included REST data proxy avoids the additional dependencies and resource setup
required by the VICC Variation Normalizer Python API. Applications can subclass
``VariationNormalizerDataProxy`` to use that API directly.
"""

import logging
import re
from abc import ABC, abstractmethod
from os import getenv

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

# Terms that are known to be unsupported by the VICC Variation Normalizer.
UNSUPPORTED_CHANGE_TERMS = frozenset(
    {
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
)


@dataclass
class MolecularProfileNameComponents:
    """Parsed molecular profile name components."""

    gene: str
    change: str | None
    c_change: str | None


class VariationNormalizerError(Exception):
    """Raised when VICC Variation Normalizer processing fails."""


class VariationNormalizerDataProxy(ABC):
    """Abstraction for VICC Variation Normalizer backends.

    Subclasses implement :meth:`normalize` to use a VICC Variation Normalizer backend.
    :meth:`normalize_molecular_profile` provides the shared CIViC molecular profile
    parsing and eligibility checks.
    """

    @abstractmethod
    def normalize(
        self,
        expr: str,
    ) -> Allele | CopyNumberChange | None:
        """Normalize a variation expression with the configured VICC backend.

        :param expr: Variation expression to normalize.
        :return: Normalized VRS variation object if normalization succeeds.
            Otherwise ``None``.
        :raises VariationNormalizerError: If the request fails or the response is not
            a valid VRS variation object.
        """
        raise NotImplementedError

    @staticmethod
    def _parse_molecular_profile_name(
        molecular_profile_name: str,
    ) -> MolecularProfileNameComponents | None:
        """Parse components from a molecular profile name.

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
    def _is_normalizable_change(
        change: str,
        molecular_profile_id: int,
    ) -> bool:
        """Determine whether a change is eligible for normalization.

        :param change: Variant change component extracted from the molecular
            profile name. This is what will be passed to the VICC Variation Normalizer.
        :param molecular_profile_id: CIViC Molecular Profile ID. Used for logging.
        :return: ``True`` if the change does not match a known unsupported pattern.
            Otherwise ``False``. This is an eligibility filter, not a guarantee that
            the VICC backend can normalize the expression.
        """
        change_lower = change.lower()

        if (
            change_lower.endswith("fs")
            or any(c in change_lower for c in ("-", "/"))
            or bool(set(change_lower.split()) & UNSUPPORTED_CHANGE_TERMS)
        ):
            _logger.warning(
                "Unsupported molecular profile change for VICC Variation Normalizer. mpid=%i, change='%s'",
                molecular_profile_id,
                change,
            )
            return False

        return True

    def _build_normalization_query(
        self,
        molecular_profile: MolecularProfile,
    ) -> str | None:
        """Build a VICC Variation Normalizer query from a molecular profile.

        :param molecular_profile: CIViC Molecular Profile.
        :return: Query to use with the VICC Variation Normalizer when the profile
            contains exactly one gene variant and a supported protein-level change.
            Otherwise ``None``.
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

        components = self._parse_molecular_profile_name(mp_name)

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

        if not self._is_normalizable_change(
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
        :raises VariationNormalizerError: If the configured VICC backend cannot
            complete the normalization request.
        """
        expression = self._build_normalization_query(molecular_profile)
        if not expression:
            return None

        normalized_variation = self.normalize(expression)

        if not normalized_variation:
            _logger.warning(
                "VICC Variation Normalizer failed to normalize query. mpid=%i, query='%s'",
                molecular_profile.id,
                expression,
            )
            return None

        return normalized_variation


class VariationNormalizerRESTDataProxy(VariationNormalizerDataProxy):
    """REST-backed data proxy for the VICC Variation Normalizer service.

    By default, requests are sent to a local service at :attr:`DEFAULT_BASE_URL`.
    Provide ``base_url`` or set :attr:`BASE_URL_ENV_VAR` to use another endpoint.
    """

    DEFAULT_BASE_URL = "http://127.0.0.1:8000/variation"
    BASE_URL_ENV_VAR = "CIVICPY_VARIATION_NORMALIZER_URL"

    def __init__(self, base_url: str | None = None):
        """Initialize the VICC Variation Normalizer REST data proxy.

        :param base_url: Base URL for the VICC Variation Normalizer REST service.
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

        :param variation: Variation object returned by the VICC Variation Normalizer.
        :param query: Original normalization query. Used for logging.
        :return: Validated VRS variation object if supported. Otherwise ``None``.
        :raises VariationNormalizerError: If the variation cannot be validated as a
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
            msg = f"VICC Variation Normalizer returned invalid VRS variation object. query={query!r}, variation_type={variation_type!r}"
            _logger.exception(msg)
            raise VariationNormalizerError(msg) from e

    def normalize(
        self,
        expr: str,
    ) -> Allele | CopyNumberChange | None:
        """Normalize a variation expression via the ``/normalize`` endpoint.

        :param expr: Variation expression to normalize.
        :return: Normalized VRS variation object if normalization succeeds.
            Otherwise ``None``.
        :raises VariationNormalizerError: If the request fails or the response is not
            a valid VRS variation object.
        """
        try:
            response = requests.get(
                f"{self.base_url}/normalize",
                params={"q": expr},
                timeout=15,
            )
            response.raise_for_status()

        except requests.exceptions.HTTPError as e:
            status_code = e.response.status_code if e.response is not None else None
            msg = (
                "VICC Variation Normalizer returned unexpected HTTP status. "
                f"query={expr!r}, "
                f"status_code={status_code}"
            )
            _logger.exception(msg)
            raise VariationNormalizerError(msg) from e

        except requests.exceptions.RequestException as e:
            msg = (
                f"VICC Variation Normalizer request failed. query={expr!r}, error={e!s}"
            )
            _logger.exception(msg)
            raise VariationNormalizerError(msg) from e

        try:
            data = response.json()
        except ValueError as e:
            msg = f"VICC Variation Normalizer returned invalid JSON. query={expr!r}"
            _logger.exception(msg)
            raise VariationNormalizerError(msg) from e

        if not isinstance(data, dict):
            msg = (
                "VICC Variation Normalizer returned an invalid response payload. "
                f"query={expr!r}"
            )
            _logger.error(msg)
            raise VariationNormalizerError(msg)

        variation = data.get("variation")
        if variation is None:
            _logger.warning(
                "VICC Variation Normalizer returned no variation object. query=%r",
                expr,
            )
            return None

        if not isinstance(variation, dict):
            msg = (
                "VICC Variation Normalizer returned an invalid variation object. "
                f"query={expr!r}"
            )
            _logger.error(msg)
            raise VariationNormalizerError(msg)

        variation = variation.copy()
        variation["name"] = expr
        variation.pop("extensions", None)

        return self._validate_vrs_variation(
            variation,
            expr,
        )
