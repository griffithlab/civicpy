"""Shared utilities for CIViC GKS exports."""

import json
from typing import Any


def serialize_canonical_json(value: Any) -> str:
    """Serialize a JSON-compatible value with stable formatting.

    :param value: Value to serialize.
    :return: Compact JSON with keys in lexical order.
    """
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
