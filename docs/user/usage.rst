.. py:module:: civicpy.civic

The **civic** module
======================

CIViCpy is primarily designed to enable exploration of the content of CIViC through Python :class:`CivicRecord` objects.
:class:`Gene`, :class:`Variant`, :class:`Assertion`, :class:`Source`, and :class:`Evidence`
are all subclasses of :class:`CivicRecord`. While these object can be generated locally, the **civic** module provides
several routines for `getting records`_.

The **civic** module may be imported from **civicpy** at the top level::

	>>>from civicpy import civic

The CivicRecord classes
-----------------------

As a base class, :class:`CivicRecord` is used to define the characteristic of all records in CIViC. This class is not
intended to be invoked directly by the end user, but provided for documentation of shared methods and variables in
child classes.

.. class:: CivicRecord

	The following methods are provided for each :class:`CivicRecord`:

	.. method:: update(allow_partial=True, force=False, **kwargs)

		Updates the record object from the cache or the server. The `allow_partial` flag will
		update the record according to the contents of CACHE, without requiring all attributes to be assigned. The
		`force` flag is used to force an update from the server, even if a full record exists in the cache. Keyword
		arguments may be passed to `kwargs`, which will update the corresponding attributes of the
		:class:`CivicRecord` instance.

	.. method:: site_link()

		Returns a URL to the record on the CIViC web application.

	Each class defines two sets of keys describing the record properties.

.. class:: Gene

.. class:: Variant

.. class:: Evidence

.. class:: Assertion

.. class:: Source

Attributes
----------

The :class:`Attribute` class is a special type of CivicRecord that is not indexed, and is used as a base container
class for additional complex records beyond those mentioned above (e.g. diseases, drugs). Attributes are not cached
except as attached objects to non-:class:`Attribute` :class:`CivicRecord` objects, and cannot be retrieved
independently.

Getting records
---------------

Records can be obtained by ID through a collection of functions provided in the `civic` module. :class:`Gene`
objects can be queried by the following methods:

.. function:: get_genes_by_ids(gene_id_list)
	A list of CIViC gene IDs are provided as `gene_id_list` and queried against the cache and (as needed) CIViC.
	Returns a list of :class:`Gene` objects.

.. function:: get_gene_by_id(gene_id)
	Similar to :func:`get_genes_by_ids`, but only one ID is passed (not a list) and only one
	:class:`Gene` returned.

.. function:: get_all_genes()
	Queries CIViC for all genes and returns as list of :class:`Gene` objects.
	The cache is not considered by this function.

.. function:: get_all_gene_ids()
	Queries CIViC for a list of all gene IDs. Useful for passing to :func:`get_genes_by_id` to
	first check cache for any previously queried genes.

Analogous methods exist for :class:`Variant`, :class:`Assertion`, :class:`Source`, and :class:`Evidence`:

.. function:: get_variants_by_ids(variant_id_list)
.. function:: get_variant_by_id(variant_id)
.. function:: get_all_variants()
.. function:: get_all_variant_ids()

.. function:: get_assertions_by_ids(assertion_id_list)
.. function:: get_assertion_by_id(assertion_id)
.. function:: get_all_assertions()
.. function:: get_all_assertion_ids()

.. function:: get_sources_by_ids(source_id_list)
.. function:: get_source_by_id(source_id)
.. function:: get_all_sources()
.. function:: get_all_source_ids()

.. function:: get_evidences_by_ids(evidence_id_list)
.. function:: get_evidence_by_id(evidence_id)
.. function:: get_all_evidences()
.. function:: get_all_evidence_ids()
