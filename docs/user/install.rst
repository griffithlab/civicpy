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

A local copy of the CIViCpy cache is kept by default at your $HOME/.civicpy/cache.pkl directory.
You may configure an environment variable CIVICPY_CACHE_FILE to adjust this destination path for
storing the cache. For example::

   >> .bashrc << echo "export CIVICPY_CACHE_FILE=$HOME/.civicpy/cache.pkl"

When loading a cache, it will be updated automatically by CIViCpy if stale (default: after 7 days).
Cache stale time is configurable when calling `civic.update_cache()`.
