.. module:: exports

The **exports** module
======================

CIViCpy supports exporting of CIViC records to Variant Call Format (VCF) files.
This enables downstream analyses such as integrating with `IGV`_, `VEP`_, and
other common bioinformatics tools. VCF exports are maintained via the :mod:`civic_vcf_writer`
and :mod:`civic_vcf_record` modules in the civicpy.exports namespace::

    >>>from civicpy.exports.civic_vcf_writer import CivicVcfWriter
    >>>from civicpy.exports.civic_vcf_record import CivicVcfRecord

CIViCpy can also export CIViC Assertions as Global Alliance for Genomics and
Health (GA4GH) Genomic Knowledge Standards (GKS) JSON. The
:mod:`civic_gks_record` module builds GKS records, while
:mod:`civic_gks_writer` writes dereferenced GKS JSON or referenced GKS Bundle
JSON::

    >>>from civicpy.exports.civic_gks_writer import CivicGksWriter
    >>>from civicpy.exports.civic_gks_record import CivicGksClinSigAssertion, CivicGksOncogenicAssertion

Other file formats are planned for future releases. Suggestions are welcome on our
`GitHub issues page <https://github.com/griffithlab/civicpy/issues>`_.

VCF
---

VCFs are written using the :class:`civicpy.exports.CivicVcfWriter` class, to
which you add :class:`civicpy.exports.CivicVcfRecord` by adding them during
initialization. :class:`civic.Variant` records can be converted to
:class:`civicpy.exports.CivicVcfRecord` records by passing in the variant
during initialization.

In order to verify whether a variant can be converted to a CivicVcfRecord object, the convenience method
``is_valid_for_vcf`` can be called on a :class:`civic.Variant` object.

Each CivicVcfRecord object passed to the CivicVcfWriter is written to the VCF file.
If two records share the same chromosome, start position, and
reference allele(s), they will not be combined into one VCF record but will instead be written
as separate VCF records. Additional CIViC data are added to the VCF as annotations to the
``CSQ`` (consequence) ``INFO`` field. All CIViC molecular profiles that the underlying variant is a part
are identified and the evidence items and assertions linked to these molecular profiles are added
to the CSQ field with one CSQ entry for each evidence item and/or assertion. Whether
a specific CSQ entry reflects an evidence item or an assertion is determined by the
``CIViC Entity Type`` CSQ field. By utilizing the CSQ field for annotations, the resulting VCF is compatible for
import into Google BigQuery (git.io/bigquery-variant-annotation).

The status of the Assertions and EvidenceItems added to the CSQ annotations can
be controlled by the ``include_status`` parameter. Only items matching the
desired include_status(es) will be added to the CSQ annotation.

.. rubric:: VCF CSQ Field Attributes
.. list-table::
   :widths: 20 70 10
   :header-rows: 1

   * - CSQ Field
     - Description
     - Compound Field [*]_
   * - Allele
     - Alternate allele
     - No
   * - Consequence
     - CIViC sequence ontology variant types for this variant
     - Yes
   * - SYMBOL
     - HGNC gene symbol for the gene associated with this variant
     - No
   * - Entrez Gene ID
     - Entrez gene identifier for the gene associated with this variant
     - No
   * - Feature_type
     - "transcript"
     - No
   * - Feature
     - The Ensembl identifier for the CIViC representative transcripts of this variant
     - No
   * - HGVSc
     - Variant representation using HGVS notation (DNA level), corresponding to the Feature
     - No
   * - HGVSp
     - Variant representation using HGVS notation (Protein level), corresponding to the Feature
     - No
   * - CIViC Variant Name
     - The CIViC variant name of this variant
     - No
   * - CIViC Variant ID
     - The CIViC internal identifier for this variant
     - No
   * - CIViC Variant Aliases
     - CIViC aliases for this variant
     - Yes
   * - CIViC Variant URL
     - CIViC URL for this variant
     - No
   * - CIViC Molecular Profile Name
     - The CIViC molecular profile name for the molecular profile of the evidence item or assertion described in this CSQ record. The molecular profile may either be a simple molecular profile for just this variant or a complex molecular profile involving this variant in combination with other CIViC variants.
     - No
   * - CIViC Molecular Profile ID
     - The CIViC internal identifier for the molecular profile
     - No
   * - CIViC Molecular Profile Aliases
     - CIViC aliases for this molecular profile
     - Yes
   * - CIViC Molecular Profile URL
     - CIViC URL for this molecular profile
     - No
   * - CIViC HGVS
     - CIViC HGVS strings for this variant
     - Yes
   * - Allele Registry ID
     - The allele registry identifier for this variant
     - No
   * - ClinVar IDs
     - ClinVar IDs associated with this variant
     - Yes
   * - CIViC Molecular Profile Score
     - The CIViC score reflecting the reelative abundance of total available curated evidence for this molecular profile
     - No
   * - CIViC Entity Type
     - The type of entity being annotated, either "evidence" or "assertion"
     - No
   * - CIViC Entity ID
     - The CIViC internal identifier for the entity being annotated
     - No
   * - CIViC Entity URL
     - The CIViC direct URL to the entity being annotated
     - No
   * - CIViC Entity Source
     - For evidence entities, the identifier of the publication used to create the evidence including the source type in the format "sourceId_(sourceType)"
     - No
   * - CIViC Entity Variant Origin
     - The variant origin of the entity being annotated, either "Somatic", "Rare Germline", "Common Germline", "Unknown", "N/A", or "Mixed"
     - No
   * - CIViC Entity Status
     - The status of the CIViC entity being annotated, either "submitted", "accepted", or "rejected"
     - No
   * - CIViC Entity Significance
     - The type of signifiance of the entity being annotated
     - No
   * - CIViC Entity Direction
     - The direction of the significance of the entity being annotated, either "Supports", or "Does Not Support"
     - No
   * - CIViC Entity Disease
     - The cancer or cancer subtype context for the entity being annotated
     - No
   * - CIViC Entity Therapies
     - A list of therapies applicable to the entity being annotated
     - Yes
   * - CIViC Entity Therapy Interaction Type
     - A term describing now more than one therapy interact with each other in the context of the entity being annotated, either "Combination", "Sequential", or "Substitutes"
     - No
   * - CIViC Evidence Phenotypes
     - A list of HPO phenotype terms linked to entity being annotated
     - Yes
   * - CIViC Evidence Level
     - For evidence entities, a level describing the robustness of the of the study supporting the evidence
     - No
   * - CIViC Evidence Rating
     - For evidence entities, a 1-5 rating indicating the curator's confidence in the quality of the summarized evidence as a number of stars
     - No
   * - CIViC Assertion ACMG Codes
     - For assertion entities, a list of ACMG codes used in the assessment of the variant under the ACMG/AMP classification guidelines
     - Yes
   * - CIViC Assertion AMP Category
     - For assertion entities, a clinical classification by AMP/ASCO/CAP guidelines
     - No
   * - CIViC Assertion NCCN Guideline
     - For assertoin entities, a string of the NCCN guideline and version
     - No
   * - CIVIC Assertion Regulatory Approval
     - For assertion entities, a boolean indicating whether or not the therapies in this assertion have regulatory approval for the treatment of the assertion disease
     - No
   * - CIVIC Assertion FDA Companion Test
     - For assertion entities, a boolean indication whether or not theassertion has an associated FDA companion test
     - No

.. [*] Compound fields contain multiple values and use the ampersand (&) character to delineate values

CivicVcfRecord
~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_vcf_record.CivicVcfRecord
   :members:

CivicVcfWriter
~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_vcf_writer.CivicVcfWriter
   :members:

Example
~~~~~~~

Here's an example of how to export all variants from CIViC to VCF::

	from civicpy import civic
    from civicpy.exports.civic_vcf_writer import CivicVcfWriter
    from civicpy.exports.civic_vcf_record import CivicVcfRecord

    records = []
    for variant in civic.get_all_variants():
        if variant.is_valid_for_vcf():
            records.append(CivicVcfRecord(variant))
    CivicVcfWriter("civic_variants.vcf", records)

.. _`IGV`: https://software.broadinstitute.org/software/igv/
.. _`VEP`: https://useast.ensembl.org/info/docs/tools/vep/index.html

GKS JSON
--------

Use :class:`civicpy.exports.civic_gks_writer.CivicGksWriter` to write
:class:`civicpy.exports.civic_gks_record.CivicGksClinSigAssertion` and
:class:`civicpy.exports.civic_gks_record.CivicGksOncogenicAssertion` records.

``CivicGksWriter`` supports two relationship representations:

* ``bundle=False`` (the default) writes **dereferenced GKS JSON**. Related GKS
  objects are included in each Assertion record.
* ``bundle=True`` writes **referenced GKS Bundle JSON**. Shared objects appear
  in keyed root collections and relationships use JSON Pointers such as
  ``#/molecularProfile/civic.mpid:33``.

.. important::

   The CIViC GKS Bundle Format is experimental in ``0.1.0``. Check
   ``bundle_format_version`` when reading a bundle.

Choose the format when creating the writer::

    # Dereferenced GKS JSON (the default)
    CivicGksWriter(Path("civic-gks-inline.json"), records)

    # Referenced GKS Bundle JSON
    CivicGksWriter(Path("civic-gks-bundle.json"), records, bundle=True)

The CLI uses separate commands because each format selects Assertions
differently::

    # One ClinVar-ready submission type as dereferenced GKS JSON
    civicpy create-gks-json --organization-id 1 --submission-type clinical_impact -o civic-gks.json

    # Both supported Statement types as referenced GKS Bundle JSON
    civicpy create-gks-bundle -o civic-gks-bundle.json

    # JSON Schema for the referenced GKS Bundle Format
    civicpy create-gks-bundle-schema

``create-gks-json`` writes ClinVar-ready Assertions approved by the specified
CIViC organization. ``--submission-type`` is required because the ClinVar
Submission API accepts one submission type at a time.

``create-gks-bundle`` writes all accepted clinical significance and
oncogenicity Assertions, even when they are not ready for ClinVar. It also
includes accepted Approvals and their organizations when available. Each
contribution stores its organization's ``clinvarAccession`` when one is
available. A Statement uses ``clinvarAccession`` for one accession and
``clinvarAccessions`` for multiple unique accessions.

``create-gks-bundle-schema`` writes the JSON Schema for the referenced bundle,
including collection key patterns and the concrete GKS models accepted by each
collection. It writes
``civic-gks-bundle-v<bundle-format-version>.schema.json`` in the current directory
and does not query CIViC or require a Variation Normalizer service.

Use ``--organization-id`` to include only Assertions approved by one
organization::

    civicpy create-gks-bundle --organization-id 1 -o civic-gks-bundle.json

What's in a GKS Bundle
~~~~~~~~~~~~~~~~~~~~~~

A bundle is one JSON object with keyed root collections. Shared variants,
genes, sources, and propositions are stored once and linked with local JSON
Pointers.

The ``metadata`` identifies the bundle format and version.
``statistics.collections`` gives the size of every collection. Collections
that can contain several object types also include counts by type.
For ``variant``, ``count`` is the number of CIViC variants, while ``types``
counts their VRS representations as Alleles or copy-number objects.

The root collections are:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Collection
     - Contents
   * - ``sequenceReference``
     - VRS sequence references, keyed by ``refgetAccession``.
   * - ``location``
     - Reusable VRS locations, including sequence locations.
   * - ``variant``
     - VRS representations grouped by ``civic.vid``, then by protein, coding,
       genomic, or unclassified coordinate level. Each representation retains
       its GA4GH identifier.
   * - ``molecularProfile``
     - CIViC molecular profiles represented as GKS categorical variants.
   * - ``feature``, ``disease``, ``phenotype``, and ``conditionSet``
     - Gene and clinical concepts used by propositions, plus grouped conditions.
   * - ``therapy`` and ``therapyGroup``
     - Individual therapies and multi-therapy groups referenced by clinical
       significance propositions.
   * - ``variantOrigin``
     - Variant origin concepts referenced by propositions. Existing external IRI
       references remain unchanged.
   * - ``source``, ``method``, and ``organization``
     - Sources and provenance objects, including the PMID-identified CIViC SOP
       and CCV framework sources.
   * - ``proposition``
     - Assertion propositions and evidence-line target propositions.
   * - ``evidence`` and ``assertion``
     - CIViC Evidence and Assertion Statements, keyed by ``civic.eid:`` and
       ``civic.aid:`` identifiers, respectively.

Objects with source IDs keep them. Variant origins use the code from their first
mapping. Propositions, condition sets, and therapy groups use a digest of their
identifying fields in the ``civic.proposition``, ``civic.conditionSet``, and
``civic.therapyGroup`` namespaces.
Descriptive changes, such as renaming a therapy or changing an alias, do not
change an enclosing group's identifier. Member and mapping order also does not
affect an identifier. Other value objects without stable identities remain
inline.

.. note::

   The GKS Bundle Format currently uses JSON Pointers in organization and
   proposition fields that VA-Spec does not type as ``iriReference``. Resolve
   these fields according to the GKS Bundle Format.

``SequenceReference`` is the only collection whose values do not carry an ``id``
field; its collection key and ``refgetAccession`` must match.

Applications can compute these IDs without creating a bundle. The function
accepts supported Pydantic GKS objects and does not modify them::

    from civicpy.exports.civic_gks_identifier import compute_civic_gks_identifier

    therapy_group_id = compute_civic_gks_identifier(therapy_group)

If records contain different representations of the same object, the bundle
keeps the first and logs a warning. Statements and collection keys are sorted so
the result is deterministic.

Each CIViC variant groups its defining VRS object and normalized members.
Molecular-profile constraints and members point to the corresponding nested VRS
representations, preserving protein, coding, and genomic HGVS forms.
A defining Allele without an HGVS expression is treated as protein-level because
CIViC uses it for the molecular profile's protein-sequence consequence.

This excerpt uses computed CIViC GKS identifiers from the test data and
shows only the Statement collections and ``proposition`` so the links are
easy to see::

    {
      "proposition": {
        "civic.proposition:-AKWXtNluL_XZYk5cDaaV7bKw6fKlPmD": {
          "id": "civic.proposition:-AKWXtNluL_XZYk5cDaaV7bKw6fKlPmD",
          "type": "VariantClinicalSignificanceProposition",
          "subjectVariant": "#/molecularProfile/civic.mpid:33",
          "geneContextQualifier": "#/feature/civic.gid:19",
          "predicate": "hasClinicalSignificanceFor",
          "objectCondition": "#/disease/civic.did:8"
        },
        "civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu": {
          "id": "civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu",
          "type": "VariantTherapeuticResponseProposition",
          "subjectVariant": "#/molecularProfile/civic.mpid:33",
          "predicate": "predictsSensitivityTo",
          "objectTherapeutic": "#/therapy/civic.tid:146"
        }
      },
      "assertion": {
        "civic.aid:6": {
          "id": "civic.aid:6",
          "type": "Statement",
          "proposition": "#/proposition/civic.proposition:-AKWXtNluL_XZYk5cDaaV7bKw6fKlPmD",
          "hasEvidenceLines": [
            {
              "type": "EvidenceLine",
              "hasEvidenceItems": ["#/evidence/civic.eid:2997"],
              "targetProposition": "#/proposition/civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu"
            }
          ]
        }
      },
      "evidence": {
        "civic.eid:2997": {
          "id": "civic.eid:2997",
          "type": "Statement",
          "proposition": "#/proposition/civic.proposition:lzu38uLu_bvAPfb7Jo_ol8741OJaSdnu"
        }
      }
    }

Assertion Statements use ``civic.aid:{id}`` keys, while Evidence Statements use
``civic.eid:{id}`` keys. Evidence lines remain inline because they have no
stable CIViC ID, but ``hasEvidenceItems`` points to supporting Statements. Follow
each Statement's ``proposition`` pointer to find its interpretation. The full
bundle links each proposition to its variant, feature, disease or phenotype, and therapy
objects in the other root collections.

To check whether an Assertion can be transformed into a
``CivicGksClinSigAssertion`` or ``CivicGksOncogenicAssertion``, call
``is_valid_for_gks_json`` on the :class:`civic.Assertion` object.

.. important::

   GKS JSON export uses the `VICC Variation Normalizer
   <https://github.com/cancervariants/variation-normalization/>`_ to build
   categorical variants. The CLI and default Python integration use its REST API,
   which requires a running service.

   The recommended setup is the `Docker installation
   <https://github.com/cancervariants/variation-normalization/#docker-installation-preferred>`_.
   By default, CIViCpy expects the service to be available at
   ``http://127.0.0.1:8000/variation``. To use a different endpoint, set the
   ``CIVICPY_VARIATION_NORMALIZER_URL`` environment variable.
   Python callers may instead pass a custom ``VariationNormalizerDataProxy`` backed
   by the VICC Variation Normalizer Python API.

CivicGksClinSigAssertion
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_gks_record.CivicGksClinSigAssertion
   :members:
   :show-inheritance:

CivicGksOncogenicAssertion
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_gks_record.CivicGksOncogenicAssertion
   :members:
   :show-inheritance:

CivicGksWriter
~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_gks_writer.CivicGksWriter
   :members:

Variation Normalization
~~~~~~~~~~~~~~~~~~~~~~~

CIViC simple molecular profiles can be normalized to GA4GH VRS Allele or Copy
Number Change objects using the `VICC Variation Normalizer`_.

By default, CIViCpy lazily creates one shared
``VariationNormalizerRESTDataProxy`` when normalization is first needed. It uses
``http://127.0.0.1:8000/variation`` by default. Set
``CIVICPY_VARIATION_NORMALIZER_URL`` before constructing the first GKS molecular
profile to use a different REST endpoint.

The backend can instead be configured to use a downstream Python implementation.
Subclass ``VariationNormalizerDataProxy``, implement ``normalize``, and configure
an instance before creating any GKS molecular profiles::

   from civicpy.exports.civic_gks_record import CivicGksMolecularProfile
   from civicpy.exports.variation_normalizer import VariationNormalizerDataProxy


   class PythonVariationNormalizer(VariationNormalizerDataProxy):
       def normalize(self, expr: str):
           return python_api_normalizer.normalize(expr)

   CivicGksMolecularProfile.configure_variation_normalizer(
       PythonVariationNormalizer()
   )

The configured backend is reused by all subsequently constructed
``CivicGksMolecularProfile`` instances. The base data proxy handles the shared
profile parsing and eligibility checks.

.. _VICC Variation Normalizer: https://github.com/cancervariants/variation-normalization/

VariationNormalizerDataProxy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: civicpy.exports.variation_normalizer.VariationNormalizerDataProxy
   :members:

VariationNormalizerRESTDataProxy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: civicpy.exports.variation_normalizer.VariationNormalizerRESTDataProxy
   :members:

Examples
~~~~~~~~

GKS records share a single variation normalizer data proxy. Set
``CIVICPY_VARIATION_NORMALIZER_URL`` before running these examples to use an endpoint
other than the default.

Here's an example of how to export all assertions to GKS JSON::

    from pathlib import Path

    from civicpy import civic
    from civicpy.exports.civic_gks_record import (
        CivicGksRecordError,
        create_gks_record_from_assertion,
    )
    from civicpy.exports.civic_gks_writer import CivicGksWriter

    records = []

    for assertion in civic.get_all_assertions():
        if assertion.is_valid_for_gks_json():
            try:
                gks_record = create_gks_record_from_assertion(
                    assertion,
                )
            except CivicGksRecordError:
                continue
            else:
                records.append(gks_record)

    CivicGksWriter(Path("gks.json"), records)

Here's an example of how to export all assertions approved by a specific organization that are
ready for submission to ClinVar.::

    from pathlib import Path

    from civicpy import civic
    from civicpy.exports.civic_gks_record import (
        CivicGksRecordError,
        create_gks_record_from_assertion,
    )
    from civicpy.exports.civic_gks_writer import CivicGksWriter

    records = []
    organization_id = 1

    for approval in civic.get_all_approvals_ready_for_clinvar_submission_for_org(organization_id):
        assertion = approval.assertion

        if assertion.is_valid_for_gks_json():
            try:
                gks_record = create_gks_record_from_assertion(
                    assertion,
                    approval=approval,
                )
            except CivicGksRecordError:
                continue
            else:
                records.append(gks_record)

    CivicGksWriter(Path("gks.json"), records)
