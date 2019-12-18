.. py:module:: civic

The **civic** module
======================

CIViCpy is primarily designed to enable exploration of the content of CIViC through Python :class:`CivicRecord` objects.
While these record objects can be initialized independently, the **civic** module also provides several routines for
`getting records`_ directly from CIViC. Use of these routines is recommended.

The **civic** module may be imported from **civicpy** at the top level::

   >>>from civicpy import civic

CIViC records
-------------

.. autoclass:: CivicRecord
   :members:

   .. automethod:: __init__

   .. attribute:: type

      The type of record. This field is set automatically by child classes and should not be changed.

   .. attribute:: id

      The record ID. This is set on initialization using the `id` keyword argument, and reflects the primary ID for
      the record as stored in CIViC.

CIViC record types
~~~~~~~~~~~~~~~~~~

The primary CIViC records are found on the CIViC advanced search page, and are fully-formed

.. autoclass:: Gene
   :show-inheritance:

   .. attribute:: description

           A curated summary of the clinical significance of this gene.

   .. attribute:: entrez_id

           The `Entrez ID`_ associated with this gene.

   .. attribute:: name

           The `HGNC Gene Symbol`_ associated with this gene.

   .. attribute:: aliases

           A list of alternate gene symbols by which this gene is referenced.

   .. attribute:: variants

           A list of :class:`Variant` records associated with this gene.

.. _Entrez ID: https://www.ncbi.nlm.nih.gov/gene/

.. _HGNC Gene Symbol: https://www.genenames.org/

.. autoclass:: Variant
   :show-inheritance:

   .. attribute:: allele_registry_id

      The `allele registry id`_ associated with this variant.

   .. attribute:: civic_actionability_score

      The CIViC `actionability score`_ associated with this variant.

   .. attribute:: description

      A curated summary of the clinical significance of this variant.

   .. attribute:: entrez_id

      The `Entrez ID`_ of the gene this variant belongs to.

   .. attribute:: entrez_name

      The `HGNC Gene Symbol`_ of the gene this variant belongs to.

   .. attribute:: gene

      The :class:`Gene` this variant belongs to.

   .. attribute:: gene_id

      The :attr:`CivicRecord.id` of the gene this variant belongs to.

   .. attribute:: name

      The curated name given to this variant.

   .. attribute:: assertions

      A list of :class:`Assertion` records associated with this variant.

   .. attribute:: clinvar_entries

      A list of `clinvar ids`_ associated with this variant.

   .. attribute:: coordinates

      A :class:`CivicAttribute` object describing `CIViC coordinates`_.

   .. attribute:: evidence_items
      evidence

      A list of :class:`Evidence` associated with this variant.

   .. attribute:: evidence_sources

      A list of sources  associated with the evidence from this variant.

   .. attribute:: hgvs_expressions

      Curated `HGVS expressions`_ describing this variant.

   .. attribute:: sources

      A list of sources objects associated with the variant description.

   .. attribute:: variant_aliases
      aliases

      A curated list of aliases by which this variant is referenced.

   .. attribute:: variant_groups
      groups

   A list of `variant groups`_ to which this variant belongs.

   .. attribute:: variant_types
      types

      A list of :class:`CivicAttribute` objects describing `variant types`_ from the
      `Sequence Ontology`_.

   .. attribute:: lifecycle_actions

      A :class:`LifecycleAction` container.

.. _allele registry id: http://reg.clinicalgenome.org

.. _actionability score: https://civicdb.org/help/variants/actionability-score

.. _clinvar ids: https://www.ncbi.nlm.nih.gov/clinvar

.. _CIViC coordinates: https://civicdb.org/help/variants/variants-coordinates

.. _HGVS expressions: https://varnomen.hgvs.org

.. _variant groups: https://civicdb.org/help/variant-groups/overview

.. _variant types: https://civicdb.org/help/variants/variants-type

.. _Sequence Ontology: http://www.sequenceontology.org/

.. autoclass:: Evidence
   :show-inheritance:

   .. attribute:: lifecycle_actions

      A :class:`LifecycleAction` container.

.. autoclass:: Assertion
   :show-inheritance:

   .. attribute:: lifecycle_actions

      A :class:`LifecycleAction` container.

CIViC attributes
~~~~~~~~~~~~~~~~

The :class:`CivicAttribute` class is a special type of CivicRecord that is not indexed, and is used as a base
class for additional complex records beyond those mentioned above (e.g. diseases, drugs). CivicAttributes are not cached
except as attached objects to non-:class:`CivicAttribute` :class:`CivicRecord` objects, and cannot be retrieved
independently.

.. autoclass:: CivicAttribute

Record provenance
~~~~~~~~~~~~~~~~~

The :class:`LifecycleAction` class is used to track the provenance of many types of :class:`CivicRecord`, by serving as
a container class for :class:`BaseLifecycleAction` objects which in turn specify the timestamp and :class:`User` associated with
a given action on a record.

.. autoclass:: LifecycleAction
   :show-inheritance:

.. autoclass:: BaseLifecycleAction
   :show-inheritance:

   .. attribute:: timestamp

      A CIViC timestamp string indicating the time the BaseLifecycleAction took place.

   .. attribute:: user

      A CIViC :class:`User` object.

.. autoclass:: Submitted
   :show-inheritance:

.. autoclass:: LastModified
   :show-inheritance:

.. autoclass:: LastReviewed
   :show-inheritance:

.. autoclass:: Accepted
   :show-inheritance:

.. autoclass:: User
   :show-inheritance:

Getting records
---------------

By ID
~~~~~

Records can be obtained by ID through a collection of functions provided in the `civic` module. :class:`Gene`
objects can be queried by the following methods:

.. autofunction:: get_gene_by_id
.. autofunction:: get_genes_by_ids
.. autofunction:: get_all_genes
.. autofunction:: get_all_gene_ids

Analogous methods exist for :class:`Variant`, :class:`Assertion`, and :class:`Evidence`:

.. autofunction:: get_variants_by_ids
.. autofunction:: get_variant_by_id
.. autofunction:: get_all_variants
.. autofunction:: get_all_variant_ids

.. autofunction:: get_assertions_by_ids
.. autofunction:: get_assertion_by_id
.. autofunction:: get_all_assertions
.. autofunction:: get_all_assertion_ids

.. autofunction:: get_all_evidence
.. autofunction:: get_all_evidence_ids

By Coordinate
~~~~~~~~~~~~~

Variant records can be searched by GRCh37 coordinates. To query specific genomic coordinates, you will
need to construct a :class:`CoordinateQuery` object, and pass this query to the
:func:`search_variants_by_coordinates` function. If you wish to query multiple genomic coordinates (e.g.
a set of variants observed in a patient tumor), construct a sorted list of :class:`CoordinateQuery` objects
(sorted by `chr`, `start`, `stop`, `alt`), and pass the list to the :func:`bulk_search_variants_by_coordinates`
function.

.. autoclass:: CoordinateQuery
.. autofunction:: search_variants_by_coordinates
.. autofunction:: bulk_search_variants_by_coordinates
