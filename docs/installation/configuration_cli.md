# Configuration

Before VariantValidator can be used, it must be configured so that it can locate the required databases and reference sequence repository.

Configuration is performed using the interactive configuration utility installed with VariantValidator.

---

# Running the configuration utility

Launch the configuration utility by running:

```bash
variantvalidator_configure
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
|---------|---------|
| `mysql` | Validator MySQL database. |
| `postgres` | VVTA PostgreSQL database. |
| `seqrepo` | Local SeqRepo installation. |
| `logging` | Logging configuration. |
| `EntrezID` | NCBI Entrez configuration. |

For example, to configure only the MySQL settings:

```bash
variantvalidator_configure --section mysql
```

or

```bash
variantvalidator_configure -s mysql
```

Similarly, to configure only the SeqRepo settings:

```bash
variantvalidator_configure --section seqrepo
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

---

# Verifying the configuration

After configuration has been completed, verify that VariantValidator can connect to the configured databases by running the test suite as described in the [Installation Guide](installation.md).

If any configuration problems are encountered, refer to the [Configuration Troubleshooting Guide](configuration_troubleshooting.md).

