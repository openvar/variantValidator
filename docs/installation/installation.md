---

# Configure MySQL for testing

The complete VariantValidator test suite creates a large number of simultaneous database connections. Before running the tests, we recommend increasing the MySQL `max_connections` setting.

Edit your MySQL configuration file (typically `my.cnf` or `mysqld.cnf`) and set:

```ini
[mysqld]
max_connections = 1000
```

Restart the MySQL server after making this change.

If you do not increase the connection limit, the test suite may fail with MySQL connection errors, particularly when running tests in parallel.

---

# Configuration

Before VariantValidator can be used, the installation must be configured.

Configuration is described in the [Configuration Guide](configuration_cli.md).

---

# Verify the installation

After completing the installation and configuration, we recommend running the VariantValidator test suite to verify that all software components and databases have been installed correctly.

If you are using the Docker databases (Quick Start installation), the test suite can be run in parallel:

```bash
pytest -n 4
```

If you performed a fully native installation, run the test suite using a single process:

```bash
pytest
```

!!! note

    Running the test suite in parallel against a native database installation may exhaust the available MySQL connections unless `max_connections` has been increased appropriately.

The test suite performs functional tests covering the VariantValidator software, the Validator database, SeqRepo, VVTA and external dependencies.

A successful installation should complete the test suite without any failures.

---

# Troubleshooting

If you encounter installation or configuration problems, see the [Configuration Troubleshooting Guide](configuration_troubleshooting.md).

---

# Next steps

Once the installation has been verified, continue to the User Manual for information on the available interfaces, supported input formats and usage examples.

