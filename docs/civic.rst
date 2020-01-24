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

   .. attribute:: lifecycle_actions

      A :class:`LifecycleAction` container.

.. _Entrez ID: https://www.ncbi.nlm.nih.gov/gene/

.. _HGNC Gene Symbol: https://www.genenames.org/

.. autoclass:: Variant
   :show-inheritance:

   .. attribute:: allele_registry_id

      The `ClinGen Allele Registry ID`_ associated with this variant.

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

.. _ClinGen Allele Registry ID: http://reg.clinicalgenome.org

.. _actionability score: https://docs.civicdb.org/en/latest/model/variants/evidence_score.html

.. _clinvar ids: https://www.ncbi.nlm.nih.gov/clinvar

.. _CIViC coordinates: https://docs.civicdb.org/en/latest/model/variants/coordinates.html

.. _HGVS expressions: https://varnomen.hgvs.org

.. _variant groups: https://docs.civicdb.org/en/latest/model/variant_groups.html

.. _variant types: https://docs.civicdb.org/en/latest/model/variants/types.html

.. _Sequence Ontology: http://www.sequenceontology.org/

.. autoclass:: Evidence
   :show-inheritance:

   .. todo:: Finish documenting

   .. attribute:: lifecycle_actions

      A :class:`LifecycleAction` container.

.. autoclass:: Assertion
   :show-inheritance:

   .. attribute:: acmg_codes

      Evidence codes used in the assessment of variants under the `ACMG/AMP`_ classification guidelines.

   .. attribute:: allele_registry_id

      The `ClinGen Allele Registry ID`_ associated with the assertion's variant.

   .. attribute:: amp_level

      The clinical interpretation classification by `AMP/ASCO/CAP`_ or `ACMG/AMP`_ guidelines.

   .. attribute:: clinical_significance

      A string indicating the type of clinical significance statement being made, values are defined based on
      the corresponding :attr:`evidence_type`. Please see `Understanding Clinical Significance`_ for more
      details on the expected values for this field.

   .. attribute:: description

      The Assertion Description gives detail including practice guidelines and approved tests for the variant.
      See `curating assertions`_ for more details.

   .. attribute:: disease

      A disease :class:`CivicAttribute`, linked to a corresponding `Disease Ontology`_ term when applicable.

   .. attribute:: drugs

      Zero or more drug :class:`CivicAttribute`, linked to corresponding `NCIT`_ terms when applicable. Only used with
      therapeutic response predictive :attr:`evidence_type`.

   .. attribute:: drug_interaction_type

      One of 'Combination', 'Sequential', or 'Substitutes', this field describes how multiple indicated drugs within
      a therapeutic response predictive assertion are related.

   .. attribute:: evidence_direction

      An indicator of whether the evidence statement supports or refutes the clinical significance of an event.

   .. attribute:: evidence_type

      Category of clinical action/relevance implicated by event. Refer to the additional `documentation on evidence types`_
      for details on how to enter evidence of each of the four types: Predictive, Prognostic, Predisposing and Diagnostic.

   .. attribute:: fda_companion_test

      A boolean indicating whether or not the assertion has an associated FDA companion test.

   .. attribute:: fda_regulatory_approval

      A boolean indicating whether or not the drugs indicated in the assertion have regulatory approval for use in
      the treatment of the assertion disease.

   .. attribute:: lifecycle_actions

      A :class:`LifecycleAction` container.

   .. attribute:: name

      A system-generated unique identifier for the assertion, e.g. `AID7`.

   .. attribute:: nccn_guideline

      A string linking the assertion to the corresponding `NCCN Guidelines for treatment of cancer by disease site`_, if
      applicable.

   .. attribute:: nccn_guideline_version

      The version associated with the indicated :attr:`nccn_guideline` document.

   .. attribute:: phenotypes

      Zero or more phenotype :class:`CivicAttribute`, linked to corresponding Human Phenotype Ontology (`HPO`_) terms
      when applicable.

   .. attribute:: status

      One of ['accepted', 'rejected', or 'submitted'], describing the state of this assertion in the CIViC curation cycle.

      - *submitted*: This assertion has been submitted by a CIViC curator or editor
      - *accepted*: This assertion has been reviewed and approved by a CIViC editor
      - *rejected*: This assertion has been reviewed and rejected by a CIViC editor

   .. attribute:: summary

      The Assertion Summary restates the Clinical Significance as a brief single sentence statement. It is intended for
      potential use in clinical reports. The Assertion Summary is designed for rapid communication of the Clinical
      Significance, especially when displayed in a longer list with other variants.

   .. attribute:: variant

      The :class:`Variant` associated with this assertion.

   .. attribute:: variant_origin

      The origin of this variant, one of ['Somatic', 'Rare Germline', 'Common Germline', 'Unknown', 'N/A', 'Germline or Somatic']

.. _AMP/ASCO/CAP: https://www.ncbi.nlm.nih.gov/pubmed/27993330

.. _ACMG/AMP: https://www.ncbi.nlm.nih.gov/pubmed/25741868

.. _curating assertions: https://docs.civicdb.org/en/latest/curating/assertions.html

.. _Disease Ontology: http://disease-ontology.org/

.. _documentation on evidence types: https://docs.civicdb.org/en/latest/model/evidence/type.html

.. _NCIT: https://ncit.nci.nih.gov/ncitbrowser/

.. _NCCN Guidelines for treatment of cancer by disease site: https://www.nccn.org/professionals/physician_gls/default.aspx#site

.. _HPO: https://hpo.jax.org/

.. _Understanding Clinical Significance: https://docs.civicdb.org/en/latest/model/evidence/clinical_significance.html#understanding-clinical-significance

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

   .. attribute:: name

      The user-defined full name for the :class:`User`.

   .. attribute:: username

      The user-defined system name for the :class:`User`.

   .. attribute:: role

      The CIViC role held by a :class:`User`: Administrator, Editor, or Curator.

   .. attribute:: area_of_expertise

      An optional attribute for a :class:`User` indicating their perspective as a CIViC participant:
         Research Scientist, Clinical Scientist, or Patient Advocate.

   .. attribute:: orcid

      An optional attribute for a :class:`User` indicating their `ORCiD <https://orcid.org/>`_.

   .. attribute:: display_name

      This attribute is populated with the first non-blank value from :attribute:`username`,
         :attribute:`name`, :attribute:`email`, or :attribute:`id`.

   .. attribute:: created_at

      A `datetime` describing when the :class:`User` object was created.

   .. attribute:: url

      A string describing a personal website or blog for the :class:`User`.

   .. attribute:: twitter_handle

      A string describing the Twitter handle for a :class:`User`.

   .. attribute:: facebook_profile

      A string describing the Facebook profile ID for a :class:`User`.

   .. attribute:: linkedin_profile

      A string describing the Linked profile ID for a :class:`User`.

   .. attribute:: bio

      A short biography for a :class:`User`.

   .. attribute:: featured_expert

      A flag indicating if a user is a featured expert, and thus displayed on the CIViC
         `domain experts <https://docs.civicdb.org/en/latest/about/domain-experts.html>`_ section

.. autoclass:: Organization
   :show-inheritance:

   .. attribute:: name

      The :class:`Organization` name.

   .. attribute:: url

      A URL for more information about the :class:`Organization`.

   .. attribute:: description

      A short text description about the :class:`Organization`.

   .. attribute:: profile_image

      A set of resource paths for the :class:`Organization` image at varying resolution.

   .. attribute:: parent

      A parent :class:`Organization`, if applicable.

.. autoclass:: Country
   :show-inheritance:

   .. attribute:: iso

      The `ISO 3166 Country Code <https://www.iso.org/iso-3166-country-codes.html>`_ for the
      :class:`Country`.

   .. attribute:: name

      The full-text name for the :class:`Country`.


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
