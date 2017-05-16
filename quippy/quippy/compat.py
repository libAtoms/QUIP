"""compat.py

Shared code for compatibility between Python versions.

"""

import sys

# Python 2.7
if sys.version_info[0] < 3:
    string_types = basestring
else:
    string_types = str

__all__ = (string_types, )
