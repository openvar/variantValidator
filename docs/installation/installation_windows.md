# Windows Installation

VariantValidator is developed and tested primarily on Linux. For Windows users, we recommend installing VariantValidator using one of the following methods.

---

# Quick Start (Recommended)

The recommended installation method is to use Docker to install the required databases while installing VariantValidator itself within Windows Subsystem for Linux (WSL2).

This approach combines the simplicity of pre-configured Docker databases with a native Linux development environment, providing full access to the Python API and command-line tools.

First install Windows Subsystem for Linux (WSL2) by following the official Microsoft installation guide:

https://learn.microsoft.com/windows/wsl/install

Once WSL2 has been installed, follow the **Quick Start** section of the [Docker Installation Guide](docker.md).

---

# Full Docker installation

The entire VariantValidator software stack, including the application and all required databases, can be run using Docker Desktop for Windows.

This approach requires minimal local configuration and is recommended for users who simply wish to deploy and run VariantValidator.

See the **Full Docker Installation** section of the [Docker Installation Guide](docker.md).

---

# Native installation under WSL2

If you prefer to install all software components yourself, VariantValidator can be installed natively within WSL2.

After installing WSL2, follow the standard Linux [Installation Guide](installation.md).

This installs VariantValidator together with the Validator, SeqRepo and VVTA databases directly within the Linux environment.

---

# Which installation method should I choose?

| Installation method | Recommended for |
|---------------------|-----------------|
| Quick Start | Most users, development, Python API and command-line tools. |
| Full Docker | Users who simply wish to deploy and run VariantValidator. |
| Native installation | Advanced users who wish to install and manage all software components themselves. |