# Common Configuration Errors

This guide describes some of the most common configuration problems encountered when installing and configuring VariantValidator, together with their typical causes and recommended solutions.

---

## Invalid SeqRepo database location

If the configured SeqRepo database location is incorrect or inaccessible, VariantValidator may fail during initialisation with errors similar to:

```text
Exception ignored in: <bound method UTA_postgresql.__del__ of <vvhgvs.dataproviders.uta.UTA_postgresql object at 0x...>>
...
AttributeError: 'UTA_postgresql' object has no attribute '_pool'
```

or

```text
NameError: name 'vval' is not defined
```

These secondary errors occur because the `Validator` object could not be created successfully.

### Cause

The `seqrepo_location` entry in the VariantValidator configuration file does not point to a valid SeqRepo installation.

### Solution

Verify that the configured SeqRepo directory exists and contains a valid SeqRepo database.

See the [Installation Guide](docker.md) for details of configuring SeqRepo.

---

## Running against MariaDB or older MySQL servers

When using MariaDB or older MySQL server versions, you may encounter an error similar to:

```text
mysql.connector.errors.NotSupportedError:
MySQL version 5.7.2 and earlier does not support COM_RESET_CONNECTION.
```

### Cause

Older MySQL-compatible servers do not support the `COM_RESET_CONNECTION` command used by newer versions of `mysql-connector-python`.

### Solution

Install version `8.0.12` of `mysql-connector-python`:

```bash
pip install "mysql-connector-python==8.0.12" --force-reinstall
```

After reinstalling the connector, verify the installation by running the VariantValidator test suite:

```bash
pytest
```

For background information, see the relevant Stack Overflow discussion:

https://stackoverflow.com/questions/58044497/is-there-a-way-to-use-pool-reset-connection-from-mysql-connector-python-with-m