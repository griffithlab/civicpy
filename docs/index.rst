.. civicpy documentation master file, created by
   sphinx-quickstart on Sat Sep  8 21:52:00 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CIViCpy
=======

v\ |release|.

**CIViCpy** is an analysis toolkit and SDK for extracting and analyzing knowledge from the `CIViC knowledgebase`_.

.. note:: The use of **Python 3** is *mandatory* for CIViCpy. This is in anticipation of the Python 2 scheduled end-of-life date \
   of `January 1, 2020 <https://www.python.org/dev/peps/pep-0373/>`_.

-------------------

**Using CIViCpy**::

   >>> from civicpy import civic
   >>> my_variant_ids = [12, 306, 79]
   >>> my_variants = civic.get_variants_by_ids(my_variant_ids)
   >>> my_variants
   [<CIViC variant 12>, <CIViC variant 79>, <CIViC variant 306>]

   >>> [(v.gene.name, v.name) for v in my_variants]
   [('BRAF', 'V600E'), ('KRAS', 'G12D'), ('ERBB2', 'AMPLIFICATION')]

   >>> braf_id = my_variants[0].gene_id # or my_variants[0].gene.id
   >>> braf_variants = civic.get_gene_by_id(braf_id).variants
   >>> len(braf_variants)
   67

   >>> set(my_variants) & set(braf_variants)
   {<CIViC variant 12>}

**CIViCpy** lets you pull data from CIViC using the :mod:`civic` module and interact with records as dynamic objects.
With the aid of caching, it is easy to explore relationships between CIViC records (e.g. assertions, genes, variants)
without worrying about querying CIViC more than once for the same data.

Features
--------

- Record Caching
- Hashable and Comparable Records
- Simplified API Requests
- Sensible Pre-caching
- Arbitrary-depth Record Nesting

The User Guide
--------------

This documentation describes how to get started with CIViCpy, the SDK, and making queries to the CIViC
knowledgebase.

.. toctree::
   :maxdepth: 2

   user/intro
   user/install
   user/civic

.. _`CIViC knowledgebase`: https://civicdb.org