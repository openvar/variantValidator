# Installation Guide

This guide describes how to install, configure and verify a VariantValidator installation.

Whether you choose a native installation or Docker deployment, these guides will help you install the required software, configure the databases and confirm that VariantValidator is operating correctly.

## Installation guides

### [Installation](installation.md)

Install VariantValidator natively on Linux or macOS using the recommended Conda environment.

This guide covers:

- Native installation
- Installing the required databases
- Verifying the installation
- Running the test suite

---

### [Docker Installation](docker.md)

Install VariantValidator using Docker.

The guide includes both:

- **Quick Start** — Docker databases with a local VariantValidator installation.
- **Full Docker Installation** — the complete VariantValidator software stack running in Docker.

---

### [Windows Installation](installation_windows.md)

Guidance for installing VariantValidator on Windows using Windows Subsystem for Linux (WSL2) or Docker Desktop.

---

## Configuration

### [Configuration Guide](configuration_cli.md)

Configure VariantValidator using the interactive configuration utility.

Topics include:

- Running the configuration utility
- Configuring database connections
- Configuring SeqRepo
- Updating existing configurations

---

### [Configuration Troubleshooting](configuration_troubleshooting.md)

Solutions to common installation and configuration problems, including database connectivity, SeqRepo configuration and dependency issues.

---

## Recommended installation path

For most users we recommend the following workflow:

1. Follow the **Quick Start** section of the Docker Installation Guide to install the required databases.
2. Install VariantValidator locally using the native Installation Guide.
3. Configure VariantValidator using the Configuration Guide.
4. Verify the installation by running the VariantValidator test suite.
5. Continue to the User Manual to begin using VariantValidator.
