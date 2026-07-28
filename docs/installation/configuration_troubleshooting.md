<img src="../static/img/logos/VV_logo.png" width="20%" />

# Common Configuration Errors

This guide describes some of the most common configuration problems encountered when installing and configuring VariantValidator, together with their typical causes and recommended solutions.

The examples below assume that VariantValidator has been installed using the VariantValidator conda environment and that the `vvenv` environment is active.

## See also

- [Installation Guide](installation.md) — Install VariantValidator and its required databases.
- [Docker Installation](docker.md) — Recommended installation using Docker-hosted databases.
- [Configuration Guide](configuration_cli.md) — Configure VariantValidator.
- [Windows Installation](installation_windows.md) — Installation and configuration under Windows.

---

# Invalid SeqRepo database location

If the configured SeqRepo database location is incorrect or inaccessible, VariantValidator may fail during initialisation with errors similar to:

```text
Exception ignored in: <bound method UTA_postgresql.__del__ of <vvhgvs.dataproviders.uta.UTA_postgresql object at 0x...>>
...
AttributeError: 'UTA_postgresql' object has no attribute '_pool'
```

or:

```text
NameError: name 'vval' is not defined
```

These secondary errors occur because the `Validator` object could not be created successfully.

## Cause

The `seqrepo_location` entry in the VariantValidator configuration file does not point to a valid SeqRepo installation.

## Solution

Verify that the configured SeqRepo directory exists and contains a valid SeqRepo database.

If you are using the Docker Quick Start installation, verify that the SeqRepo data have been extracted to the location specified in the [Docker Installation Guide](docker.md).

If you performed a native installation, see the SeqRepo section of the [Installation Guide](installation.md).

After correcting the SeqRepo location, ensure that the `vvenv` conda environment is active and restart VariantValidator.

---

# Running against MariaDB or older MySQL servers

When using MariaDB or older MySQL server versions, you may encounter an error similar to:

```text
mysql.connector.errors.NotSupportedError:
MySQL version 5.7.2 and earlier does not support COM_RESET_CONNECTION.
```

## Cause

Older MySQL-compatible servers do not support the `COM_RESET_CONNECTION` command used by newer versions of `mysql-connector-python`.

## Solution

Ensure that the VariantValidator conda environment is active:

```bash
conda activate vvenv
```

Install version `8.0.12` of `mysql-connector-python`:

```bash
pip install "mysql-connector-python==8.0.12" --force-reinstall
```

After reinstalling the connector, verify the installation by running the VariantValidator test suite:

```bash
pytest
```

For additional background information, see the [Stack Overflow discussion of `pool_reset_connection` and mysql-connector-python](https://stackoverflow.com/questions/58044497/is-there-a-way-to-use-pool-reset-connection-from-mysql-connector-python-with-m).

---

# Getting help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter installation or configuration problems that are not resolved by the guidance above, we encourage you to seek assistance.

Before contacting the development team, you may find the following documentation helpful:

- [Installation Guide](installation.md)
- [Docker Installation](docker.md)
- [Configuration Guide](configuration_cli.md)
- [Windows Installation](installation_windows.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>

