.. py:module:: civic

The **civic** module
======================

CIViCpy is primarily designed to enable exploration of the content of CIViC through Python :class:`CivicRecord` objects.
While these record objects can be initialized independently, the **civic** module also provides several routines for
:ref:`getting_records` directly from CIViC. Use of these routines is recommended.

The **civic** module may be imported from **civicpy** at the top level::

   >>>from civicpy import civic

CIViC Records
-------------

.. autoclass:: CivicRecord
   :members:

   .. automethod:: __init__

   .. attribute:: type

      The type of record. This field is set automatically by child classes and should not be changed.

   .. attribute:: id

      The record ID. This is set on initialization using the `id` keyword argument, and reflects the primary ID for
      the record as stored in CIViC.

The primary CIViC records are found on the sidebar menu on CIViC, and are fully-formed.

Gene
^^^^

.. autoclass:: Gene
   :show-inheritance:
   :members:

   .. attribute:: aliases

      A list of alternate gene symbols by which this gene is referenced.

   .. attribute:: description

      A curated summary of the clinical significance of this gene.

   .. attribute:: entrez_id

      The `Entrez ID`_ associated with this gene.

   .. attribute:: name

      The `HGNC Gene Symbol`_ associated with this gene.

.. _Entrez ID: https://www.ncbi.nlm.nih.gov/gene/

.. _HGNC Gene Symbol: https://www.genenames.org/

Factor
^^^^^^

.. autoclass:: Factor
   :show-inheritance:
   :members:

   .. attribute:: aliases

      A list of alternate names by which this factor is referenced.

   .. attribute:: description

      A curated summary of the clinical significance of this factor.

   .. attribute:: full_name

      Factor names are often an commonly-used abbreviation. The full name is
      the unabbreviated name.

   .. attribute:: name

      The shortest, most concise reference to the factor. Often an
      abbreviation.

   .. attribute:: ncit_id

      The `NCIthesaurus ID`_ referencing the factor.

.. _NCIthesaurus ID: https://ncithesaurus.nci.nih.gov/ncitbrowser/

Fusion
^^^^^^

.. autoclass:: Fusion
   :show-inheritance:
   :members:

   .. attribute:: subtype

   .. attribute:: aliases

      A list of alternate names by which this fusion is referenced.

   .. attribute:: description

      A curated summary of the clinical significance of this fusion.

   .. attribute:: five_prime_gene_id

      The :attr:`CivicRecord.id` of the 5' fusion partner :class:`Gene` if that partner is
      ``KNOWN``.

   .. attribute:: five_prime_partner_status

      The status of the 5' fusion partner. One of ``KNOWN``, ``UNKNOWN``, or
      ``MULTIPLE``.

   .. attribute:: name

      The name of the fusion. This will be the 5' partner, followed by the 3'
      partner, separated by ``::``. If a partner is ``KNOWN``, the `HGNC Gene Symbol`_
      of the partner gene is used. If the partner is ``UNKNOWN``,
      a ``?`` is used. If there are ``MULTIPLE`` possible gene partners,
      ``v`` is used.

   .. attribute:: three_prime_gene_id

      The :attr:`CivicRecord.id` of the 3' fusion partner :class:`Gene` if that partner is
      ``KNOWN``.

   .. attribute:: three_prime_partner_status

      The status of the 3' fusion partner. One of ``KNOWN``, ``UNKNOWN``, or
      ``MULTIPLE``.

.. _HGNC Gene Symbol: https://www.genenames.org/


Variant
^^^^^^^

.. autoclass:: Variant
   :show-inheritance:
   :members:

   .. attribute:: feature_id

      The :attr:`CivicRecord.id` of the :class:`Gene`, :class:`Factor`, or
      :class:`Fusion` the variant belongs to.

   .. attribute:: name

      The curated name given to this variant.

   .. attribute:: single_variant_molecular_profile_id

      The :attr:`CivicRecord.id` of the :class:`MolecularProfile` representing the single
      variant on its own.

   .. attribute:: subtype

      The specific type of variant. One of ``gene_variant``,
      ``factor_variant``, or ``fusion_variant``.

   .. attribute:: variant_aliases

      A curated list of aliases by which this variant is referenced.

   .. attribute:: variant_types

      A list of :class:`CivicAttribute` objects describing `variant types`_ from the
      `Sequence Ontology`_.

.. _variant types: https://docs.civicdb.org/en/latest/model/variants/types.html

.. _Sequence Ontology: http://www.sequenceontology.org/


GeneVariant
"""""""""""

.. autoclass:: GeneVariant
   :show-inheritance:
   :members:

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

   .. attribute:: hgvs_expressions

      Curated `HGVS expressions`_ describing this variant.

.. _ClinGen Allele Registry ID: http://reg.clinicalgenome.org

.. _clinvar ids: https://www.ncbi.nlm.nih.gov/clinvar

.. _CIViC coordinates: https://docs.civicdb.org/en/latest/model/variants/coordinates.html

.. _HGVS expressions: https://varnomen.hgvs.org


FactorVariant
"""""""""""""

.. autoclass:: FactorVariant
   :show-inheritance:
   :members:

   .. attribute:: ncit_id

      The `NCIthesaurus ID`_ referencing the factor variant.

.. _NCIthesaurus ID: https://ncithesaurus.nci.nih.gov/ncitbrowser/


FusionVariant
"""""""""""""

.. autoclass:: FusionVariant
   :show-inheritance:
   :members:

   .. attribute:: five_prime_coordinates

      A :class:`CivicAttribute` object describing `CIViC coordinates`_ of the
      5' fusion partner, if that partner is ``KNOWN``.

   .. attribute:: three_prime_coordinates

      A :class:`CivicAttribute` object describing `CIViC coordinates`_ of the
      3' fusion partner, if that partner is ``KNOWN``.

   .. attribute:: vicc_compliant_name

      A name representing the fusion variant compliant with the `VICC fusion
      specification`_.

.. _CIViC coordinates: https://docs.civicdb.org/en/latest/model/variants/coordinates.html

.. _VICC fusion specification: https://fusions.cancervariants.org/en/latest/nomenclature.html


MolecularProfile
^^^^^^^^^^^^^^^^

.. autoclass:: MolecularProfile
   :show-inheritance:
   :members:

   .. attribute:: aliases

      A curated list of aliases by which this molecular profile is referenced.

   .. attribute:: description

      A curated summary of the clinical significance of this molecular
      profile.

   .. attribute:: molecular_profile_score

      The CIViC `molecular profile score`_ associated with this molecular
      profile.

   .. attribute:: name

      The human readable name of this molecular profile, including gene and variant names.

   .. attribute:: source_ids

      A list of integers designating the :attr:`CivicRecord.id` for the
      class:`Source` records associated with the molecular profile
      description.

   .. attribute:: variant_ids

      An list of integers designating the :attr:`CivicRecord.id` for the class:`Variant` records involved in this
      molecular profile.

.. _molecular profile score: https://civic.readthedocs.io/en/latest/model/molecular_profiles/evidence_score.html


Evidence
^^^^^^^^

.. autoclass:: Evidence
   :show-inheritance:
   :members:

   .. attribute:: assertion_ids

      The list of :attr:`CivicRecord.id` of :class:`Assertion` records this evidence is a part of.

   .. attribute:: description

      The Evidence Statement (returned as `description` by the CIViC API) is a brief summary of the clinical implications of the :attr:`variant` in the context of the specific :attr:`disease`, :attr:`evidence_type`, and :attr:`significance` as curated from the cited literature source.

   .. attribute:: disease_id

      The :attr:`CivicRecord.id` of the :class:`Disease` record of the cancer of cancer subtype context for the evidence record. **None** for functional evidence_type.

   .. attribute:: evidence_direction

      One of 'Supports', 'Does Not Support', indicating whether the evidence statement supports or refutes the significance of an event.

   .. attribute:: evidence_level

      The evidence level describes the robustness of the study supporting the evidence item. Five different evidence levels are supported: “A - Validated association”, “B - Clinical evidence”, “C - Case study”, “D - Preclinical evidence”, and “E - Inferential association”. For more information, please see `Understanding Levels`_.

   .. attribute:: evidence_type

      Category of clinical action/relevance implicated by event. Refer to the additional `documentation on evidence types`_
      for details on how to enter evidence of each of the six types: Predictive, Prognostic, Predisposing, Diagnostic, Functional, and Oncogenic.

   .. attribute:: molecular_profile_id

      The :attr:`CivicRecord.id` of the :class:`MolecularProfile` this evidence item belongs to.

   .. attribute:: name

      A system-generated unique identifier for the evidence record, e.g. `EID7`.

   .. attribute:: phenotype_ids

      The list of :attr:`CivicRecord.id` of :class:`Phenotype` records linked to corresponding `Human Phenotype Ontology (HPO)`_ terms when applicable.

   .. attribute:: rating

      The Evidence Rating is an integer from 1 to 5, indicating the curator’s confidence in the quality of the summarized evidence as a number of stars. For more information about this metric, please see `Understanding Evidence Ratings`_.

   .. attribute:: significance

      A string indicating the type of significance statement being made, values are defined based on
      the corresponding :attr:`evidence_type`. Please see `Understanding Significance`_ for more
      details on the expected values for this field.

   .. attribute:: source_id

      The :attr:`CivicRecord.id` of the :class:`Source` object this evidence was derived from.

   .. attribute:: status

      One of 'accepted', 'rejected', or 'submitted', describing the state of this evidence in the CIViC curation cycle. An evidence item needs to be reviewed by a CIViC editor before being accepted or rejected. Therefore "submitted" evidence might not be accurate or complete.

      - *submitted*: This evidence has been submitted by a CIViC curator or editor
      - *accepted*: This evidence has been reviewed and approved by a CIViC editor
      - *rejected*: This evidence has been reviewed and rejected by a CIViC editor

   .. attribute:: therapy_ids

      The list of :attr:`CivicRecord.id` of the :class:`Therapy` objects this evidence item is linked to. Only used with therapeutic response predictive evidence_type.

   .. attribute:: therapy_interaction_type

      One of 'Combination', 'Sequential', or 'Substitutes', this field describes how multiple indicated therapies within
      a therapeutic response predictive :attr:`evidence_type` are related.

.. _Human Phenotype Ontology (HPO): https://hpo.jax.org/

.. _Understanding Levels: https://civic.readthedocs.io/en/latest/model/evidence/level.html#understanding-levels

.. _Understanding Evidence Ratings: https://civic.readthedocs.io/en/latest/model/evidence/evidence_rating.html#understanding-evidence-ratings


Assertion
^^^^^^^^^

.. autoclass:: Assertion
   :show-inheritance:
   :members:

   .. attribute:: acmg_codes

      Evidence codes used in the assessment of germline variant pathogenicity under the `ACMG/AMP`_ classification guidelines.

   .. attribute:: amp_level

      The clinical tiering of somatic variants by `AMP/ASCO/CAP`_ guidelines.

   .. attribute:: assertion_direction

      One of 'Supports' or 'Does Not Support', indicating whether the evidence statement supports or refutes the significance of an event.

   .. attribute:: assertion_type

      Category of clinical action/relevance implicated by event. Refer to the additional `documentation on assertion types`_
      for details on how to enter assertions of each of the five types: Predictive, Prognostic, Predisposing, Diagnostic, and Oncogenic.

   .. attribute:: clingen_codes

      Classification of somatic variant oncogenicity under the `ClinGen/CGC/VICC`_ classification guidelines.

   .. attribute:: description

      The Assertion Description gives detail including practice guidelines and approved tests for the molecular profile.
      See `curating assertions`_ for more details.

   .. attribute:: disease_id

      The :attr:`CivicRecord.id` of the :class:`Disease` record of the cancer of cancer subtype context for the assertion record.

   .. attribute:: evidence_ids

      A list of :attr:`CivicRecord.id` of the :class:`Evidence` records supporting this assertion record.

   .. attribute:: fda_companion_test

      A boolean indicating whether or not the assertion has an associated FDA companion test.

   .. attribute:: fda_regulatory_approval

      A boolean indicating whether or not the therapies indicated in the assertion have regulatory approval for use in
      the treatment of the assertion disease.

   .. attribute:: molecular_profile_id

      The :attr:`CivicRecord.id` of the molecular profile this assertion belongs to.

   .. attribute:: name

      A system-generated unique identifier for the assertion, e.g. `AID7`.

   .. attribute:: nccn_guideline

      A string linking the assertion to the corresponding `NCCN Guidelines for treatment of cancer by disease site`_, if
      applicable.

   .. attribute:: nccn_guideline_version

      The version associated with the indicated :attr:`nccn_guideline` document.

   .. attribute:: phenotype_ids

      Zero or more :class:`Phenotype` :attr:`CivicRecord.id`, linked to corresponding Human Phenotype Ontology (`HPO`_) terms
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

   .. attribute:: therapy_interaction_type

      One of 'Combination', 'Sequential', or 'Substitutes', this field describes how multiple indicated therapies within
      a therapeutic response predictive :attr:`evidence_type` assertion are related.

   .. attribute:: variant_origin

      The origin of the variants in this molecular profile, one of 'Somatic', 'Rare Germline', 'Common Germline', 'Unknown', 'N/A', 'Germline or Somatic', or 'Mixed'

.. _AMP/ASCO/CAP: https://www.ncbi.nlm.nih.gov/pubmed/27993330

.. _ACMG/AMP: https://www.ncbi.nlm.nih.gov/pubmed/25741868

.. _ClinGen/CGC/VICC: https://pubmed.ncbi.nlm.nih.gov/35101336/

.. _curating assertions: https://docs.civicdb.org/en/latest/curating/assertions.html

.. _Disease Ontology: http://disease-ontology.org/

.. _documentation on evidence types: https://docs.civicdb.org/en/latest/model/evidence/type.html

.. _documentation on assertion types: https://docs.civicdb.org/en/latest/model/assertions/overview.html

.. _NCIT: https://ncit.nci.nih.gov/ncitbrowser/

.. _NCCN Guidelines for treatment of cancer by disease site: https://www.nccn.org/professionals/physician_gls/default.aspx#site

.. _HPO: https://hpo.jax.org/

.. _Understanding Significance: https://civic.readthedocs.io/en/latest/model/evidence/significance.html#understanding-significance


Source
^^^^^^

.. autoclass:: Source
   :members:

   .. attribute:: abstract

      The abstract text of the source.

   .. attribute:: asco_abstract_id

      For ASCO sources, the abstract ID.

   .. attribute:: author_string

      A string of all of the authors of the source or, for ASCO sources, the abstract presenter.

   .. attribute:: citation

      A short string containing key information about the source for human-readable
      identification.

   .. attribute:: citation_id

      A unique identifier for the source. For PubMed sources, this is the
      PMID. For ASH sources this is the DOI. For ASCO sources this is the
      ASCO Web ID found in the URL of the abstract.

   .. attribute:: clinical_trials

      A list of `Clinical Trial`_ IDs described in the source.

   .. attribute:: full_journal_title

      The full title of the publishing journal.

   .. attribute:: journal

      An abbreviated version of the title of the publishing journal.

   .. attribute:: pmc_id

      When available, the `PubMed Central`_ ID of the source.

   .. attribute:: publication_date

      The date the source was published.

   .. attribute:: source_type

      The platform making the source available. One of ``PUBMED``, ``ASCO``,
      or ``ASH``.

   .. attribute:: source_url

      A link to the source on the platfrom that made the source available.

   .. attribute:: title

      The title of the source.

.. _Clinical Trial: https://clinicaltrials.gov/

.. _PubMed Central: https://pmc.ncbi.nlm.nih.gov/


Disease
^^^^^^^

.. autoclass:: Disease
   :members:

   .. attribute:: aliases

      A list of alternate names for the disease.

   .. attribute:: disease_url

      A link to the `Disease Ontology`_ entry for the disease concept.

   .. attribute:: doid

      The `Disease Ontology`_ ID for the disease concept.

   .. attribute:: name

      The name of the disease.

.. _Disease Ontology: http://disease-ontology.org/


Therapy
^^^^^^^

.. autoclass:: Therapy
   :members:

   .. attribute:: aliases

      A list of alternate names for the therapy.

   .. attribute:: name

      The name of the therapy.

   .. attribute:: ncit_id

      The `NCIthesaurus`_ ID for the therapy concept.

   .. attribute:: therapy_url

      A link to the `NCIthesaurus`_ entry for the therapy concept.

.. _NCIthesaurus: https://ncithesaurus.nci.nih.gov/ncitbrowser/


Phenotype
^^^^^^^^^

.. autoclass:: Phenotype
   :members:

   .. attribute:: name

      The name of the phenotype.

   .. attribute:: hpo_id

      The `Human Phenotype Ontology`_ ID for the phenotype concept.

   .. attribute:: phenotype_url

      A link to the `Human Phenotype Ontology`_ entry for the phenotype concept.

.. _Human Phenotype Ontology: https://hpo.jax.org/


CIViC Attributes
^^^^^^^^^^^^^^^^

The :class:`CivicAttribute` class is a special type of CivicRecord that is not indexed, and is used as a base
class for additional complex records beyond those mentioned above (e.g. diseases, therapies). CivicAttributes are not cached
except as attached objects to non-:class:`CivicAttribute` :class:`CivicRecord` objects, and cannot be retrieved
independently.

.. autoclass:: CivicAttribute

