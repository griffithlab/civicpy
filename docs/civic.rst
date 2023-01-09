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

   .. attribute:: site_link

      A URL string to the appropriate landing page for the CivicRecord on the CIViC web application.

CIViC record types
~~~~~~~~~~~~~~~~~~

The primary CIViC records are found on the CIViC advanced search page, and are fully-formed

.. autoclass:: Gene
   :show-inheritance:

   .. attribute:: aliases

      A list of alternate gene symbols by which this gene is referenced.

   .. attribute: assertions

      A list of :class:`Assertion` records that this gene is involved in.

   .. attribute:: description

      A curated summary of the clinical significance of this gene.

   .. attribute:: entrez_id

      The `Entrez ID`_ associated with this gene.

   .. attribute:: name

      The `HGNC Gene Symbol`_ associated with this gene.

   .. attribute:: sources

      A list of :class:`CivicAttribute` source objects associated with the gene description.

   .. attribute:: variants

      A list of :class:`Variant` records associated with this gene.

.. _Entrez ID: https://www.ncbi.nlm.nih.gov/gene/

.. _HGNC Gene Symbol: https://www.genenames.org/

.. autoclass:: Variant
   :show-inheritance:

   .. attribute:: allele_registry_id

      The `ClinGen Allele Registry ID`_ associated with this variant.

   .. attribute:: clinvar_entries

      A list of `clinvar ids`_ associated with this variant.

   .. attribute:: coordinates

      A :class:`CivicAttribute` object describing `CIViC coordinates`_.

   .. attribute:: entrez_id

      The `Entrez ID`_ of the gene this variant belongs to.

   .. attribute:: entrez_name

      The `HGNC Gene Symbol`_ of the gene this variant belongs to.

   .. attribute:: gene

      The :class:`Gene` this variant belongs to.

   .. attribute:: gene_id

      The :attr:`CivicRecord.id` of the gene this variant belongs to.

   .. attribute:: hgvs_expressions

      Curated `HGVS expressions`_ describing this variant.

   .. attribute:: name

      The curated name given to this variant.

   .. attribute:: moleulcar_profiles

      A list of :class:`MolecularProfile` objects of all the molecular
      profiles involving this variant.

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

.. _ClinGen Allele Registry ID: http://reg.clinicalgenome.org

.. _clinvar ids: https://www.ncbi.nlm.nih.gov/clinvar

.. _CIViC coordinates: https://docs.civicdb.org/en/latest/model/variants/coordinates.html

.. _HGVS expressions: https://varnomen.hgvs.org

.. _variant groups: https://docs.civicdb.org/en/latest/model/variant_groups.html

.. _variant types: https://docs.civicdb.org/en/latest/model/variants/types.html

.. _Sequence Ontology: http://www.sequenceontology.org/

.. autoclass:: MolecularProfile
   :show-inheritance:

   .. attribute: aliases

      A curated list of aliases by which this molecular profile is referenced.

   .. attribute: assertions

      A list of :class:`Assertion` records associated with this molecular
      profile.

   .. attribute::  definition

      A curated summary of the clinical significance of this molecular
      profile.

   .. attribute: evidence_items

      A list of :class:`Evidence` associated with this molecular profile.

   .. attribute:: molecular_profile_score

      The CIViC `molcular profile score`_ associated with this molecular
      profile.

   .. attribute:: name

      The human readable name of this molecular profile, including gene and variant names.

   .. attribute: variant_ids

      An list of integers designating the :attr:`CivicRecord.id` for the variants involved in this
      molecular profile.

   .. attribute: sources

      A list of :class:`CivicAttribute` source objects associated with the molecular profile description.

   .. attribute: variants

      A list :class:`Variant` objects involved in this molecular profile.

.. _molecular profile score: https://civic.readthedocs.io/en/latest/model/molecular_profiles/evidence_score.html

.. autoclass:: Evidence
   :show-inheritance:

   .. attribute:: assertions

      CIViC :class:`Assertion` records containing this evidence.

   .. attribute:: description
      statement

      The Evidence Statement (returned as `description` by the CIViC API) is a brief summary of the clinical implications of the :attr:`variant` in the context of the specific :attr:`disease`, :attr:`evidence_type`, and :attr:`significance` as curated from the cited literature source.

   .. attribute:: disease

      The cancer or cancer subtype context for the evidence record.

   .. attribute:: evidence_direction

      One of 'Supports', 'Does Not Support', indicating whether the evidence statement supports or refutes the significance of an event.

   .. attribute:: evidence_level

      The evidence level describes the robustness of the study supporting the evidence item. Five different evidence levels are supported: “A - Validated association”, “B - Clinical evidence”, “C - Case study”, “D - Preclinical evidence”, and “E - Inferential association”. For more information, please see `Understanding Levels`_.

   .. attribute:: evidence_type

      Category of clinical action/relevance implicated by event. Refer to the additional `documentation on evidence types`_
      for details on how to enter evidence of each of the six types: Predictive, Prognostic, Predisposing, Diagnostic, Functional, and Oncogenic.

   .. attribute:: molecular_profile

      The :class:`MolecularProfile` object this evidence item belongs to.

   .. attribute:: molecular_profile_id

      The :attr:`CivicRecord.id` of the molecular profile this evidence item belongs to.

   .. attribute:: name

      A system-generated unique identifier for the evidence record, e.g. `EID7`.

   .. attribute:: phenotypes

      Zero or more phenotype :class:`CivicAttribute`, linked to corresponding Human Phenotype Ontology (`HPO`_) terms
      when applicable.

   .. attribute:: rating

      The Evidence Rating is an integer from 1 to 5, indicating the curator’s confidence in the quality of the summarized evidence as a number of stars. For more information about this metric, please see `Understanding Evidence Ratings`_.

   .. attribute:: significance

      A string indicating the type of significance statement being made, values are defined based on
      the corresponding :attr:`evidence_type`. Please see `Understanding Significance`_ for more
      details on the expected values for this field.

   .. attribute:: source

      A :class:`CivicAttribute` source object from which this evidence was derived.

   .. attribute:: status

      One of 'accepted', 'rejected', or 'submitted', describing the state of this evidence in the CIViC curation cycle. An evidence item needs to be reviewed by a CIViC editor before being accepted or rejected. Therefore "submitted" evidence might not be accurate or complete.

      - *submitted*: This evidence has been submitted by a CIViC curator or editor
      - *accepted*: This evidence has been reviewed and approved by a CIViC editor
      - *rejected*: This evidence has been reviewed and rejected by a CIViC editor

   .. attribute:: therapies

      Zero or more therapy :class:`CivicAttribute`, linked to corresponding `NCIT`_ terms when applicable. Only used with
      therapeutic response predictive :attr:`evidence_type`.

   .. attribute:: therapy_interaction_type

      One of 'Combination', 'Sequential', or 'Substitutes', this field describes how multiple indicated therapies within
      a therapeutic response predictive :attr:`evidence_type` are related.

.. _Understanding Levels: https://civic.readthedocs.io/en/latest/model/evidence/level.html#understanding-levels

.. _Understanding Evidence Ratings: https://civic.readthedocs.io/en/latest/model/evidence/evidence_rating.html#understanding-evidence-ratings

.. autoclass:: Assertion
   :show-inheritance:

   .. attribute:: acmg_codes

      Evidence codes used in the assessment of variants under the `ACMG/AMP`_ classification guidelines.

   .. attribute:: amp_level

      The clinical interpretation classification by `AMP/ASCO/CAP`_ or `ACMG/AMP`_ guidelines.

   .. attribute:: assertion_direction

      One of 'Supports' or 'Does Not Support', indicating whether the evidence statement supports or refutes the significance of an event.

   .. attribute:: assertion_type

      Category of clinical action/relevance implicated by event. Refer to the additional `documentation on assertion types`_
      for details on how to enter assertions of each of the five types: Predictive, Prognostic, Predisposing, Diagnostic, and Oncogenic.

   .. attribute:: description

      The Assertion Description gives detail including practice guidelines and approved tests for the molecular profile.
      See `curating assertions`_ for more details.

   .. attribute:: disease

      A disease :class:`CivicAttribute`, linked to a corresponding `Disease Ontology`_ term when applicable.

   .. attribute:: fda_companion_test

      A boolean indicating whether or not the assertion has an associated FDA companion test.

   .. attribute:: fda_regulatory_approval

      A boolean indicating whether or not the therapies indicated in the assertion have regulatory approval for use in
      the treatment of the assertion disease.

   .. attribute:: molecular_profile

      The :class:`MolecularProfile` object this assertion belongs to.

   .. attribute:: molecular_profile_id

      The :attr:`CivicRecord.id` of the molecular profile this assertion belongs to.

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

   .. attribute:: significance

      A string indicating the type of significance statement being made, values are defined based on
      the corresponding :attr:`evidence_type`. Please see `Understanding Significance`_ for more
      details on the expected values for this field.

   .. attribute:: status

      One of 'accepted', 'rejected', or 'submitted', describing the state of this assertion in the CIViC curation cycle. An Assertion needs to be reviewed by a CIViC editor before being accepted or rejected. Therefore "submitted" Assertions might not be accurate or complete.

      - *submitted*: This assertion has been submitted by a CIViC curator or editor
      - *accepted*: This assertion has been reviewed and approved by a CIViC editor
      - *rejected*: This assertion has been reviewed and rejected by a CIViC editor

   .. attribute:: summary

      The Assertion Summary restates the Significance as a brief single sentence statement. It is intended for
      potential use in clinical reports. The Assertion Summary is designed for rapid communication of the 
      Significance, especially when displayed in a longer list with other molecular profiles.

   .. attribute:: therapies

      Zero or more therapy :class:`CivicAttribute`, linked to corresponding `NCIT`_ terms when applicable. Only used with
      therapeutic response predictive :attr:`evidence_type`.

   .. attribute:: therapy_interaction_type

      One of 'Combination', 'Sequential', or 'Substitutes', this field describes how multiple indicated therapies within
      a therapeutic response predictive :attr:`evidence_type` assertion are related.

   .. attribute:: variant_origin

      The origin of the variants in this molecular profile, one of 'Somatic', 'Rare Germline', 'Common Germline', 'Unknown', 'N/A', 'Germline or Somatic', or 'Mixed'

.. _AMP/ASCO/CAP: https://www.ncbi.nlm.nih.gov/pubmed/27993330

.. _ACMG/AMP: https://www.ncbi.nlm.nih.gov/pubmed/25741868

.. _curating assertions: https://docs.civicdb.org/en/latest/curating/assertions.html

.. _Disease Ontology: http://disease-ontology.org/

.. _documentation on evidence types: https://docs.civicdb.org/en/latest/model/evidence/type.html

.. _documentation of assertion types: https://docs.civicdb.org/en/latest/model/assertions/overview.html

.. _NCIT: https://ncit.nci.nih.gov/ncitbrowser/

.. _NCCN Guidelines for treatment of cancer by disease site: https://www.nccn.org/professionals/physician_gls/default.aspx#site

.. _HPO: https://hpo.jax.org/

.. _Understanding Significance: https://civic.readthedocs.io/en/latest/model/evidence/significance.html#understanding-significance

CIViC attributes
~~~~~~~~~~~~~~~~

The :class:`CivicAttribute` class is a special type of CivicRecord that is not indexed, and is used as a base
class for additional complex records beyond those mentioned above (e.g. diseases, therapies). CivicAttributes are not cached
except as attached objects to non-:class:`CivicAttribute` :class:`CivicRecord` objects, and cannot be retrieved
independently.

.. autoclass:: CivicAttribute

Getting records
---------------

By ID
~~~~~

Records can be obtained by ID through a collection of functions provided in the `civic` module. :class:`Gene`
objects can be queried by the following methods:

.. autofunction:: get_gene_by_id
.. autofunction:: get_genes_by_ids
.. autofunction:: get_all_genes

Analogous methods exist for :class:`Variant`, :class:`MolecularProfile`, :class:`Assertion`, and :class:`Evidence`:

.. autofunction:: get_variant_by_id
.. autofunction:: get_variants_by_ids
.. autofunction:: get_all_variants

.. autofunction:: get_molecular_profile_by_id
.. autofunction:: get_molecular_profiles_by_ids
.. autofunction:: get_all_molecular_profiles

.. autofunction:: get_assertion_by_id
.. autofunction:: get_assertions_by_ids
.. autofunction:: get_all_assertions

.. autofunction:: get_evidence_by_id
.. autofunction:: get_evidence_by_ids
.. autofunction:: get_all_evidence

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

Coordinates can also be used to query :class:`Assertion` records:

.. autofunction:: search_assertions_by_coordinates

By Other Attribute
~~~~~~~~~~~~~~~~~~

.. autofunction:: search_variants_by_allele_registry_id
.. autofunction:: search_variants_by_hgvs
.. autofunction:: search_variants_by_name
