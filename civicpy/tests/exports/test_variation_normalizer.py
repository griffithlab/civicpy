"""Test that Variation Normalizer data proxy works correctly"""

from unittest.mock import Mock, patch

import pytest
from requests.exceptions import HTTPError, RequestException

from civicpy.civic import get_molecular_profile_by_id
from civicpy.exports.variation_normalizer import (
    VariationNormalizerError,
    VariationNormalizerRestDataProxy,
)


@pytest.fixture
def rest_dp():
    return VariationNormalizerRestDataProxy()


BRAF_V600E_QUERY = "BRAF V600E"


class TestVariationNormalizerRestDataProxy(object):
    """Test that VariationNormalizerRestDataProxy works as expected

    These are mocked for CI purposes
    """

    @patch("requests.get")
    @pytest.mark.parametrize(
        ("query", "expected_vrs_fixture"),
        [
            (BRAF_V600E_QUERY, "braf_v600e_vrs"),
            ("BRAF Amplification", "braf_amplification_vrs"),
        ],
    )
    def test_normalize_success(
        self, mock_get, request, rest_dp, query, expected_vrs_fixture
    ):

        expected_vrs = request.getfixturevalue(expected_vrs_fixture)

        mock_response = Mock()
        mock_response.raise_for_status.return_value = None
        mock_response.json.return_value = {
            "variation": expected_vrs.model_dump(exclude_none=True)
        }
        mock_get.return_value = mock_response

        result = rest_dp.normalize(query)
        assert result == expected_vrs

    @patch("requests.get")
    def test_normalize_returns_none(self, mock_get, rest_dp):
        mock_response = Mock()
        mock_response.raise_for_status.return_value = None
        mock_response.json.return_value = {"variation": None}
        mock_get.return_value = mock_response

        result = rest_dp.normalize(BRAF_V600E_QUERY)
        assert result is None

    @patch("requests.get")
    def test_normalize_returns_no_variation_key(self, mock_get, caplog, rest_dp):
        mock_response = Mock()
        mock_response.raise_for_status.return_value = None
        mock_response.json.return_value = {"variant": {None}}
        mock_get.return_value = mock_response

        result = rest_dp.normalize(BRAF_V600E_QUERY)
        assert result is None
        assert (
            "Variation Normalizer returned no variation object. query='BRAF V600E'"
            in caplog.text
        )

    @patch("requests.get")
    def test_normalize_http_error(self, mock_get, rest_dp):
        mock_response = Mock()
        mock_response.raise_for_status.side_effect = HTTPError(
            response=Mock(status_code=500)
        )
        mock_get.return_value = mock_response

        with pytest.raises(
            VariationNormalizerError,
            match="Variation Normalizer returned unexpected HTTP status. query='BRAF V600E', status_code=500",
        ):
            rest_dp.normalize(BRAF_V600E_QUERY)

    @patch("requests.get")
    def test_normalize_request_error(self, mock_get, rest_dp):
        mock_get.side_effect = RequestException("Connection failed")

        with pytest.raises(
            VariationNormalizerError,
            match="Variation Normalizer request failed. query='BRAF V600E'",
        ):
            rest_dp.normalize(BRAF_V600E_QUERY)

    @patch("requests.get")
    def test_normalize_validation_error(self, mock_get, rest_dp):
        mock_response = Mock()
        mock_response.raise_for_status.return_value = None
        mock_response.json.return_value = {"variation": {"type": "Allele"}}
        mock_get.return_value = mock_response

        with pytest.raises(
            VariationNormalizerError,
            match="Variation Normalizer returned invalid VRS variation object. query='BRAF V600E', variation_type='Allele'",
        ):
            rest_dp.normalize(BRAF_V600E_QUERY)

    @patch("requests.get")
    def test_normalize_molecular_profile(
        self, mock_get, rest_dp, braf_v600e_vrs, v600e_mp
    ):
        mock_response = Mock()
        mock_response.raise_for_status.return_value = None
        mock_response.json.return_value = {
            "variation": braf_v600e_vrs.model_dump(exclude_none=True)
        }
        mock_get.return_value = mock_response

        result = rest_dp.normalize_molecular_profile(v600e_mp)
        assert result == braf_v600e_vrs

    @pytest.mark.parametrize(
        ("mp_id", "log_msg"),
        [
            (
                1615,
                "Molecular profiles containing cDNA changes are not supported. mpid=1615, name='VHL R167Q (c.500G>A)'",
            ),
            (
                332,
                "Unsupported molecular profile change for Variation Normalizer. mpid=332, change='Mutation'",
            ),
            (
                2995,
                "Variant type 'FusionVariant' is not supported. mpid=2995, name='RUNX1::RUNX1T1 Fusion'",
            ),
            (
                5117,
                "Variant type 'FactorVariant' is not supported. mpid=5117, name='TMB High'",
            ),
            (
                5879,
                "Variant type 'RegionVariant' is not supported. mpid=5879, name='1q21.2 Amplification'",
            ),
            (
                4344,
                "Complex molecular profiles are not supported. mpid=4344, name='EZH2 Y646S OR EZH2 Y646F OR EZH2 Y646H OR EZH2 Y646C OR EZH2 Y646N OR EZH2 A692V OR EZH2 A682G'",
            ),
        ],
    )
    def test_not_supported(self, caplog, rest_dp, mp_id, log_msg):
        rest_dp.normalize = Mock()

        mp = get_molecular_profile_by_id(mp_id)
        result = rest_dp.normalize_molecular_profile(mp)
        rest_dp.normalize.assert_not_called()
        assert result is None
        assert log_msg in caplog.text

    @pytest.mark.parametrize(
        "mp_id",
        [
            302,
        ],
    )
    def test_normalize_called(self, rest_dp, mp_id):
        rest_dp.normalize = Mock(return_value=None)

        mp = get_molecular_profile_by_id(mp_id)
        rest_dp.normalize_molecular_profile(mp)
        rest_dp.normalize.assert_called_once_with("ERBB2 Amplification")


@pytest.mark.live
class TestVariationNormalizerRestDataProxyLive(object):
    """Test that VariationNormalizerRestDataProxy works as expected

    Note: This hits the live Variation Normalizer service
    """

    @pytest.mark.parametrize(
        ("query", "expected_vrs_fixture"),
        [
            (BRAF_V600E_QUERY, "braf_v600e_vrs"),
            ("BRAF Amplification", "braf_amplification_vrs"),
        ],
    )
    def test_normalize(self, request, query, expected_vrs_fixture, rest_dp):
        result = rest_dp.normalize(query)
        assert result == request.getfixturevalue(expected_vrs_fixture)
