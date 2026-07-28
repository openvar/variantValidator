<img src="../static/img/logos/VV_logo.png" width="20%" />

# Windows Installation

VariantValidator is developed and tested primarily on Linux. Windows users can run VariantValidator using Windows Subsystem for Linux (WSL2) or Docker Desktop.

For most users and developers, we recommend the **Quick Start** installation. This uses WSL2 to provide the Linux environment and Docker to provide the required VariantValidator databases.

Users who prefer to run the complete VariantValidator software stack in Docker can use the **Full Docker installation**.

Advanced users who wish to install and manage all database components themselves can perform a **Native installation under WSL2**.

## See also

- [Installation Guide](installation.md) — Main installation instructions for VariantValidator.
- [Docker Installation Guide](docker.md) — Quick Start and Full Docker installation instructions.
- [Configuration Guide](configuration_cli.md) — Configure VariantValidator after installation.
- [Configuration Troubleshooting Guide](configuration_troubleshooting.md) — Resolve common installation and configuration problems.

---

# Quick Start (Recommended)

The recommended installation method for Windows is to use Docker to provide the required Validator, SeqRepo and VVTA databases while installing VariantValidator and VariantFormatter within Windows Subsystem for Linux (WSL2).

This approach combines the simplicity of the pre-configured Docker databases with a Linux development environment, providing full access to the VariantValidator and VariantFormatter Python APIs, command-line tools and development environment.

First install Windows Subsystem for Linux (WSL2) by following the [official Microsoft WSL installation guide](https://learn.microsoft.com/windows/wsl/install).

Once WSL2 has been installed, follow the **Quick Start** section of the [Docker Installation Guide](docker.md).

The VariantValidator conda environment should be created and activated within WSL2 as described in the installation instructions.

---

# Full Docker installation

The entire VariantValidator software stack, including the application and all required databases, can be run using Docker Desktop for Windows.

This approach requires minimal local configuration and is recommended for users who wish to deploy and run VariantValidator without maintaining a local development installation.

See the **Full Docker Installation** section of the [Docker Installation Guide](docker.md).

---

# Native installation under WSL2

If you prefer to install all software components yourself, VariantValidator can be installed natively within WSL2.

After installing WSL2, follow the standard Linux [Installation Guide](installation.md).

The installation uses the VariantValidator conda environment defined by `environment.yml` and installs VariantValidator and VariantFormatter together with the Validator, SeqRepo and VVTA databases within the Linux environment.

---

# Which installation method should I choose?

| Installation method | Recommended for |
|---------------------|-----------------|
| Quick Start | Most users, development, Python API and command-line tools. |
| Full Docker | Users who wish to deploy and run VariantValidator without maintaining a local development installation. |
| Native installation | Advanced users who wish to install and manage all software components themselves. |

---

# Getting help

VariantValidator has been developed to support a wide range of users, from those new to HGVS nomenclature to experienced clinical scientists and bioinformaticians. If you encounter difficulties installing or configuring VariantValidator on Windows, we encourage you to seek assistance.

Before contacting the development team, you may find the following documentation helpful:

- [Installation Guide](installation.md)
- [Docker Installation Guide](docker.md)
- [Configuration Guide](configuration_cli.md)
- [Configuration Troubleshooting Guide](configuration_troubleshooting.md)

If you still require assistance, you can contact the VariantValidator team using our [contact form](https://variantvalidator.org/help/contact/).

Software bugs and feature requests can be reported through the [VariantValidator GitHub issue tracker](https://github.com/openvar/VariantValidator/issues).

---

## Acknowledgements

**VariantValidator was originally developed at the University of Leicester (2016–2019). It is now maintained and developed by the University of Manchester, with continued hosting and development contributions from the University of Leicester.**

<img src="../static/img/logos/Manchester_logo.png" width="40%" align="left"/>
<img src="../static/img/logos/uniofleicesterlogo.png" width="40%" align="right" />
<br clear="both"/>
