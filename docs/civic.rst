The **civic** module
======================

CIViCpy is primarily designed to enable exploration of the content of CIViC through Python :class:`CivicRecord` objects.
While these record objects can be initialized independently, the **civic** module also provides several routines for
`getting records`_ directly from CIViC. Use of these routines is recommended.

The **civic** module may be imported from **civicpy** at the top level::

	>>>from civicpy import civic

CIViC records
-------------

.. autoclass:: civic.CivicRecord
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

.. autoclass:: civic.Gene
   :members:

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

.. autoclass:: civic.Variant

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

.. _allele registry id: http://reg.clinicalgenome.org

.. _actionability score: https://civicdb.org/help/variants/actionability-score

.. _clinvar ids: https://www.ncbi.nlm.nih.gov/clinvar

.. _CIViC coordinates: https://civicdb.org/help/variants/variants-coordinates

.. _HGVS expressions: https://varnomen.hgvs.org

.. _variant groups: https://civicdb.org/help/variant-groups/overview

.. _variant types: https://civicdb.org/help/variants/variants-type

.. _Sequence Ontology: http://www.sequenceontology.org/

.. autoclass:: civic.Evidence

.. autoclass:: civic.Assertion

CIViC attributes
~~~~~~~~~~~~~~~~

The :class:`CivicAttribute` class is a special type of CivicRecord that is not indexed, and is used as a base container
class for additional complex records beyond those mentioned above (e.g. diseases, drugs). CivicAttributes are not cached
except as attached objects to non-:class:`CivicAttribute` :class:`CivicRecord` objects, and cannot be retrieved
independently.

.. class:: CivicAttribute

Getting records
---------------

Records can be obtained by ID through a collection of functions provided in the `civic` module. :class:`Gene`
objects can be queried by the following methods:

.. automodule:: civic
   :members: get_gene_by_id, get_genes_by_ids, get_all_genes, get_all_gene_ids

Analogous methods exist for :class:`Variant`, :class:`Assertion`, and :class:`Evidence`:

.. automodule:: civic
   :members: get_variants_by_ids, get_variant_by_id, get_all_variants, get_all_variant_ids
   :undoc-members:
   :noindex:

.. automodule:: civic
   :members: get_assertions_by_ids, get_assertion_by_id, get_all_assertions, get_all_assertion_ids
   :undoc-members:
   :noindex:

.. automodule:: civic
   :members: get_all_evidence, get_all_evidence_ids
   :undoc-members:
   :noindex:

