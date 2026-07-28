<img src="../static/img/logos/VV_logo.png" width="20%" />

# Configuration

Before VariantValidator can be used, it must be configured so that it can locate the required databases and reference sequence repository.

Configuration is performed using the interactive configuration utility installed with VariantValidator.

## See also

- [Installation Guide](installation.md) — Install VariantValidator and its required databases.
- [Docker Installation](docker.md) — Recommended installation using Docker-hosted databases.
- [Configuration Troubleshooting](configuration_troubleshooting.md) — Resolve common configuration problems.
- [Windows Installation](installation_windows.md) — Installation and configuration under Windows.

---

# Running the configuration utility

Ensure that the VariantValidator conda environment is active:

```bash
conda activate vvenv
```

Launch the configuration utility by running:

```bash
variantvalidator-configure
```

If no user configuration file exists, a new configuration file is created automatically using the default settings supplied with VariantValidator.

The utility then steps through each configuration section, prompting for new values.

For example:

```text
Section: mysql
user [vvadmin]:
password [password]:
host [127.0.0.1]:
port [3306]:
database [validator]:
```

Press **Enter** to keep the current value, or type a replacement value.

When configuration is complete, the updated configuration is written automatically.

---

# Configuration file location

The user configuration file is stored in the VariantValidator configuration directory.

The location can be determined programmatically using:

```python
from VariantValidator import settings

print(settings.get_config_dir())
```

If no configuration file is present, the configuration utility creates one automatically from the packaged default configuration.

---

# Configuring individual sections

Rather than editing the entire configuration, individual sections may be configured independently.

The available sections are:

| Section | Purpose |
| --- | --- |
| `mysql` | Validator MySQL database. |
| `postgres` | VVTA PostgreSQL database. |
| `seqrepo` | Local SeqRepo installation. |
| `logging` | Logging configuration. |
| `EntrezID` | NCBI Entrez configuration. |

For example, to configure only the MySQL settings:

```bash
variantvalidator-configure --section mysql
```

or:

```bash
variantvalidator-configure -s mysql
```

Similarly, to configure only the SeqRepo settings:

```bash
variantvalidator-configure --section seqrepo
```

---

# Configuration sections

## MySQL

This section configures the Validator annotation database.

Typical settings include:

- Username
- Password
- Host
- Port
- Database name

These values should match the MySQL database created during installation.

---

## PostgreSQL

This section configures the VVTA transcript alignment database.

Typical settings include:

- Username
- Password
- Host
- Port
- Database name

These values should match the PostgreSQL database created during installation.

---

## SeqRepo

This section specifies the location of the local SeqRepo installation.

The configured directory should contain the extracted SeqRepo data downloaded during installation.

---

## Logging

The logging section controls the behaviour of VariantValidator logging.

For most users the default settings are appropriate and do not require modification.

---

## EntrezID

This section configures access to the NCBI Entrez services.

Most users can leave the default values unchanged.

---

# Reconfiguring VariantValidator

The configuration utility may be run as many times as required.

Existing values are displayed as defaults, allowing individual settings to be updated without modifying the remainder of the configuration.

Ensure that the VariantValidator conda environment is active before running the utility:

```bash
conda activate vvenv
```

Then run:

```bash
variantvalidator-configure
```

---

# Verifying the configuration

After configuration has been completed, verify that VariantValidator can connect to the configured databases by running the test suite as described in the [Installation Guide](installation.md).

If configuration problems are encountered, refer to the [Configuration Troubleshooting Guide](configuration_troubleshooting.md).

---

# Getting help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter difficulties configuring VariantValidator, we encourage you to seek assistance.

Before contacting the development team, you may find the following documentation helpful:

- [Installation Guide](installation.md)
- [Docker Installation](docker.md)
- [Configuration Troubleshooting](configuration_troubleshooting.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>
