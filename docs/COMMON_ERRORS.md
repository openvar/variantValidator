# Here is a file of common errors which we have seen

## Configuration

#### Invalid SeqRepo database location
An invalid SeqRepo location in your config will generate an error like this 
```python
>>> validate = vval.validate(variant, genome_build, select_transcripts)
Exception ignored in: <bound method UTA_postgresql.__del__ of <vvhgvs.dataproviders.uta.UTA_postgresql object at 0x7f4bb03bf1d0>>
Traceback (most recent call last):
  File "/local/miniconda3/envs/vvweb/lib/python3.6/site-packages/vvhgvs/dataproviders/uta.py", line 477, in __del__
    self.close()
  File "/local/miniconda3/envs/vvweb/lib/python3.6/site-packages/vvhgvs/dataproviders/uta.py", line 481, in close
    self._pool.closeall()
AttributeError: 'UTA_postgresql' object has no attribute '_pool'
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'vval' is not defined
>>> validation = validate.format_as_dict(with_meta=True)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'validate' is not defined
>>> print(json.dumps(validation, sort_keys=True, indent=4, separators=(',', ': ')))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'validation' is not defined
```
***Update your seqrepo location in your config file. See INSTALLATION.md***

#### Running on a system with MariaDB
see [StackOverflow](https://stackoverflow.com/questions/58044497/is-there-a-way-to-use-pool-reset-connection-from-mysql-connector-python-with-m)
```python
>>> validation = validate.format_as_dict(with_meta=True)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'validate' is not defined
>>> validate = vval.validate(variant, genome_build, select_transcripts)
ERROR: <class 'mysql.connector.errors.NotSupportedError'> Reset session is not supported with MySQL server version 5.7.2 or earlier.
CRITICAL: <class 'mysql.connector.errors.NotSupportedError'> Reset session is not supported with MySQL server version 5.7.2 or earlier.
Traceback (most recent call last):
  File "/local/miniconda3/envs/vvweb/lib/python3.6/site-packages/mysql/connector/connection_cext.py", line 768, in reset_session
    self.cmd_reset_connection()
  File "/local/miniconda3/envs/vvweb/lib/python3.6/site-packages/mysql/connector/connection_cext.py", line 678, in cmd_reset_connection
    raise errors.NotSupportedError("MySQL version 5.7.2 and "
mysql.connector.errors.NotSupportedError: MySQL version 5.7.2 and earlier does not support COM_RESET_CONNECTION.
```
***Solution***
Downgrade mysql-connector-python to version 8.0.12
```bash
$ pip install 'mysql-connector-python==8.0.12' --force-reinstall
```
Test by running pytest
