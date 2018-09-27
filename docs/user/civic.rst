.. py:module:: civic

The **civic** module
======================

CIViCpy is primarily designed to enable exploration of the content of CIViC through Python :class:`CivicRecord` objects.
While these record objects can be initialized independently, the **civic** module also provides several routines for
`getting records`_ directly from CIViC. Use of these routines is recommended.

The **civic** module may be imported from **civicpy** at the top level::

	>>>from civicpy import civic

The CivicRecord class
---------------------

As a base class, :class:`CivicRecord` is used to define the characteristic of all records in CIViC. This class is not
intended to be invoked directly by the end user, but provided for documentation of shared methods and variables in
child classes.

.. class:: CivicRecord

	.. method:: __init__(partial=False, **kwargs):
		The record object may be initialized by the user, though the practice is discouraged. To do so, values for each
		of the object attributes (except `type`) must be specified as keyword arguments, or the `partial` parameter must
		be set to **True**. If `partial` is set to **True**, the `id` keyword argument is still required.

		Users are encouraged to use the below functions for `getting records`_ in lieu of directly initializing record
		objects.

	.. method:: update(allow_partial=True, force=False, **kwargs)

		Updates the record object from the cache or the server. The `allow_partial` flag will
		update the record according to the contents of CACHE, without requiring all attributes to be assigned. The
		`force` flag is used to force an update from the server, even if a full record exists in the cache. Keyword
		arguments may be passed to `kwargs`, which will update the corresponding attributes of the
		:class:`CivicRecord` instance.

	.. method:: site_link()

		Returns a URL to the record on the CIViC web application.

	.. attribute:: type

		The type of record. This field is set automatically by child classes and should not be changed.

	.. attribute:: id

		The record ID. Must be set on initialization using the `id` keyword argument.

Primary CIViC Records
~~~~~~~~~~~~~~~~~~~~~

The primary CIViC records are those found on

.. class:: Gene(CivicRecord)

	.. attribute:: description

		A free-text description summarizing the known association between the gene and cancer.

	.. attribute:: entrez_id

		The `Entrez ID <https://www.ncbi.nlm.nih.gov/gene/>`_ associated with this gene.

.. class:: Variant(CivicRecord)

.. class:: Evidence(CivicRecord)

.. class:: Assertion(CivicRecord)

.. class:: Source(CivicRecord)

The CivicAttributes Class
~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`CivicAttribute` class is a special type of CivicRecord that is not indexed, and is used as a base container
class for additional complex records beyond those mentioned above (e.g. diseases, drugs). CivicAttributes are not cached
except as attached objects to non-:class:`CivicAttribute` :class:`CivicRecord` objects, and cannot be retrieved
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
