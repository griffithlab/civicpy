Installation
============

Installation is easy as `*PyPI*` with the `pip` package manager.

pip
---
To install with pip::

   >> pip install civicpy

That's it!

.. _config-cache:

Configuring Cache Save
----------------------

If you wish to persist a local copy of the CIViCpy cache, you may set an environment variable
CIVICPY_CACHE_FILE to a destination path for storing the cache. For example::

   >> .bashrc << echo "export CIVICPY_CACHE_FILE=$HOME/.civicpy/cache.pkl"

When loading a cache, it will be updated automatically by CIViCpy if stale (default: after 7 days).
Cache stale time is configurable when calling `civic.update_cache()`.
