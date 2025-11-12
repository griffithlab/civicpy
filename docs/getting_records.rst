.. py:module:: civic
   :noindex:

.. _getting_records:

Getting Records
===============

CIViCpy offeres a wide range of convenience methods as part of the `civic` module to retrieve different
CIViC entities.

Get All Records For A Specific Entity Type
------------------------------------------

Features
~~~~~~~~

.. autofunction:: get_all_features
.. autofunction:: get_all_genes
.. autofunction:: get_all_factors
.. autofunction:: get_all_fusions

Variants
~~~~~~~~

.. autofunction:: get_all_variants
.. autofunction:: get_all_gene_variants
.. autofunction:: get_all_factor_variants
.. autofunction:: get_all_fusion_variants

Molecular Profiles
~~~~~~~~~~~~~~~~~~

.. autofunction:: get_all_molecular_profiles

Assertions
~~~~~~~~~~

.. autofunction:: get_all_assertions
.. autofunction:: get_all_assertions_ready_for_clinvar_submission_for_org

Evidence Items
~~~~~~~~~~~~~~

.. autofunction:: get_all_evidence

Variant Groups
~~~~~~~~~~~~~~

.. autofunction:: get_all_variant_groups

Sources
~~~~~~~

.. autofunction:: get_all_sources

Diseases
~~~~~~~~

.. autofunction:: get_all_diseases

Therapies
~~~~~~~~~

.. autofunction:: get_all_therapies

Phenotypes
~~~~~~~~~~

.. autofunction:: get_all_phenotypes

Organizations
~~~~~~~~~~~~~

.. autofunction:: get_all_organizations

Approvals

.. autofunction:: get_all_approvals
.. autofunction:: get_all_approvals_ready_for_clinvar_submission_for_org

By ID
-----

Records can be obtained by CIViC ID through a collection of functions provided in the `civic` module.

Features
~~~~~~~~

.. autofunction:: get_feature_by_id
.. autofunction:: get_features_by_ids

.. autofunction:: get_gene_by_id
.. autofunction:: get_genes_by_ids

.. autofunction:: get_factor_by_id
.. autofunction:: get_factors_by_ids

.. autofunction:: get_fusion_by_id
.. autofunction:: get_fusions_by_ids

Variants
~~~~~~~~

.. autofunction:: get_variant_by_id
.. autofunction:: get_variants_by_ids

Molecular Profiles
~~~~~~~~~~~~~~~~~~

.. autofunction:: get_molecular_profile_by_id
.. autofunction:: get_molecular_profiles_by_ids

Assertions
~~~~~~~~~~

.. autofunction:: get_assertion_by_id
.. autofunction:: get_assertions_by_ids

Evidence Items
~~~~~~~~~~~~~~

.. autofunction:: get_evidence_by_id
.. autofunction:: get_evidence_by_ids

Variant Groups
~~~~~~~~~~~~~~

.. autofunction:: get_variant_group_by_id
.. autofunction:: get_variant_groups_by_ids

Sources
~~~~~~~

.. autofunction:: get_source_by_id
.. autofunction:: get_sources_by_ids

Diseases
~~~~~~~~

.. autofunction:: get_disease_by_id
.. autofunction:: get_diseases_by_ids

Therapies
~~~~~~~~~

.. autofunction:: get_therapy_by_id
.. autofunction:: get_therapies_by_ids

Phenotypes
~~~~~~~~~~

.. autofunction:: get_phenotype_by_id
.. autofunction:: get_phenotypes_by_ids

Organizations
~~~~~~~~~~~~~

.. autofunction:: get_organization_by_id
.. autofunction:: get_organizations_by_ids

Approvals
~~~~~~~~~

.. autofunction:: get_approval_by_id
.. autofunction:: get_approvals_by_ids


By Coordinates
--------------

Variant records can be searched by GRCh37 coordinates. To query specific genomic coordinates, you will
need to construct a :class:`CoordinateQuery` object, and pass this query to the
:func:`search_variants_by_coordinates` function. If you wish to query multiple genomic coordinates (e.g.
a set of variants observed in a patient tumor), construct a sorted list of :class:`CoordinateQuery` objects
(sorted by `chr`, `start`, `stop`, `alt`), and pass the list to the :func:`bulk_search_variants_by_coordinates`
function.

.. autoclass:: CoordinateQuery
.. autofunction:: search_variants_by_coordinates
.. autofunction:: bulk_search_variants_by_coordinates

Coordinates can also be used to query :class:`Assertion` and
:class:`Evidence` records:

.. autofunction:: search_assertions_by_coordinates
.. autofunction:: search_evidence_by_coordinates

By Other Attribute
------------------

Genes
~~~~~

.. autofunction:: get_gene_by_entrez_id
.. autofunction:: get_gene_by_name

Factors
~~~~~~~

.. autofunction:: get_factor_by_ncit_id
.. autofunction:: get_factor_by_name

Fusions
~~~~~~~

.. autofunction:: get_fusion_by_name
.. autofunction:: search_fusions_by_partner_gene_id

Variants
~~~~~~~~

.. autofunction:: search_variants_by_allele_registry_id
.. autofunction:: search_variants_by_hgvs
.. autofunction:: search_variants_by_name

Sources
~~~~~~~

.. autofunction:: get_pubmed_source_by_id
.. autofunction:: get_ash_source_by_doi
.. autofunction:: get_asco_source_by_id

Diseases
~~~~~~~~

.. autofunction:: get_disease_by_doid
.. autofunction:: get_disease_by_name

Therapies
~~~~~~~~~

.. autofunction:: get_therapy_by_ncit_id
.. autofunction:: get_therapy_by_name

Phenotypes
~~~~~~~~~~

.. autofunction:: get_phenotype_by_hpo_id
.. autofunction:: get_phenotype_by_name

Approvals
~~~~~~~~~~~~

.. autofunction:: search_approvals_by_organization_id
.. autofunction:: search_approvals_by_assertion_id
