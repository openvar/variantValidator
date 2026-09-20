import importlib
import os

import VariantValidator as VV
from VariantValidator import settings


def test_validator_database_pool_sizes(monkeypatch):
    """
    Verify that Validator uses the configured MySQL and VVTA PostgreSQL
    connection pool sizes.
    """

    # ------------------------------------------------------------------
    # Default configuration
    # ------------------------------------------------------------------
    monkeypatch.delenv(
        "VALIDATOR_MYSQL_POOL_SIZE",
        raising=False,
    )
    monkeypatch.delenv(
        "VVTA_POSTGRES_POOL_SIZE",
        raising=False,
    )

    importlib.reload(settings)

    validator = VV.Validator()

    assert validator.db.pool.pool_size == 1
    assert validator.hdp._pool.maxconn == 1

    del validator

    # ------------------------------------------------------------------
    # Environment overrides
    # ------------------------------------------------------------------
    monkeypatch.setenv(
        "VALIDATOR_MYSQL_POOL_SIZE",
        "5",
    )
    monkeypatch.setenv(
        "VVTA_POSTGRES_POOL_SIZE",
        "5",
    )

    importlib.reload(settings)

    validator = VV.Validator()

    assert validator.db.pool.pool_size == 5
    assert validator.hdp._pool.maxconn == 5

    del validator

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
