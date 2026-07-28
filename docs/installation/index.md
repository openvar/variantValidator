<img src="../static/img/logos/VV_logo.png" width="20%" />

# Installation Guide

This guide describes how to install, configure and verify a VariantValidator installation.

For most users and developers, we recommend the **Quick Start** installation, which uses Docker to provide the required databases while VariantValidator and VariantFormatter are installed locally in the VariantValidator conda environment.

Native and Full Docker installation options are also available.

## See also

- [Installation](installation.md) — Install VariantValidator and VariantFormatter on Linux or macOS.
- [Docker Installation](docker.md) — Use Docker for the VariantValidator databases or complete software stack.
- [Windows Installation](installation_windows.md) — Install using WSL2 or Docker Desktop.
- [Configuration Guide](configuration_cli.md) — Configure VariantValidator after installation.

---

# Installation guides

## Installation

[Installation](installation.md)

Install VariantValidator and VariantFormatter on Linux or macOS using the VariantValidator conda environment defined by `environment.yml`.

This guide covers:

- creating the VariantValidator conda environment;
- installing VariantValidator and VariantFormatter;
- installing the required databases;
- configuring the installation;
- verifying the installation; and
- running the test suite.

---

## Docker Installation

[Docker Installation](docker.md)

Install VariantValidator using Docker.

The guide includes two installation options:

- **Quick Start** — use Docker to provide the Validator, SeqRepo and VVTA databases while installing VariantValidator and VariantFormatter locally.
- **Full Docker Installation** — run the complete VariantValidator software stack using Docker.

The Quick Start installation is recommended for most users and developers.

---

## Windows Installation

[Windows Installation](installation_windows.md)

Install VariantValidator on Windows using Windows Subsystem for Linux (WSL2) or Docker Desktop.

The Windows guide covers the Quick Start, Full Docker and native WSL2 installation options.

---

# Configuration

## Configuration Guide

[Configuration Guide](configuration_cli.md)

Configure VariantValidator using the interactive configuration utility.

Topics include:

- running the configuration utility;
- configuring database connections;
- configuring SeqRepo; and
- updating an existing configuration.

---

## Configuration Troubleshooting

[Configuration Troubleshooting](configuration_troubleshooting.md)

Solutions to common installation and configuration problems, including database connectivity, SeqRepo configuration and dependency issues.

---

# Recommended installation path

For most users and developers, we recommend the following workflow:

1. Follow the **Quick Start** section of the [Docker Installation Guide](docker.md) to install the required databases.
2. Create the VariantValidator conda environment defined by `environment.yml`.
3. Install VariantValidator and VariantFormatter.
4. Configure VariantValidator using the [Configuration Guide](configuration_cli.md).
5. Verify the installation by running the VariantValidator test suite.
6. Continue to the User Manual to begin using VariantValidator.

---

# Getting help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter difficulties installing or configuring VariantValidator, we encourage you to seek assistance.

Before contacting the development team, you may find the following documentation helpful:

- [Installation](installation.md)
- [Docker Installation](docker.md)
- [Windows Installation](installation_windows.md)
- [Configuration Guide](configuration_cli.md)
- [Configuration Troubleshooting](configuration_troubleshooting.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>

