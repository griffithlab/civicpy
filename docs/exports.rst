.. module:: exports

The **exports** module
======================

CIViCpy supports exporting of CIViC records to Variant Call Format (VCF) files.
This enables downstream analyses such as integrating with `IGV`_, `VEP`_, and
other common bioinformatics tools. VCF exports are maintained via the :mod:`civic_vcf_writer`
and :mod:`civic_vcf_record` modules in the civicpy.exports namespace::

    >>>from civicpy.exports.civic_vcf_writer import CivicVcfWriter
    >>>from civicpy.exports.civic_vcf_record import CivicVcfRecord

CIViCPy also supports exporting of CIViC assertion records to JSON files, where
assertions are represented as Global Alliance for Genomics and Health (GA4GH) Genomic
Knowledge Standards (GKS) objects. GKS JSON exports are maintained via the
:mod:`civic_gks_writer` and :mod:`civic_gks_record` modules in the civicpy.exports
namespace::

    >>>from civicpy.exports.civic_gks_writer import CivicGksWriter
    >>>from civicpy.exports.civic_gks_record import CivicGksPredictiveAssertion, CivicGksDiagnosticAssertion, CivicGksPrognosticAssertion

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

GKS JSON files are written using the :class:`civicpy.exports.civic_gks_writer.CivicGksWriter`
class to which you add :class:`civicpy.exports.civic_gks_record.CivicGksPredictiveAssertion`,
:class:`civicpy.exports.civic_gks_record.CivicGksDiagnosticAssertion`, or
:class:`civicpy.exports.civic_gks_record.CivicGksPrognosticAssertion` during initialization.

In order to verify whether an assertion can be converted to a CivicGksAssertionRecord
object, the convenience method ``is_valid_for_gks_json`` can be called on a
:class:`civic.Assertion` object.

_CivicGksAssertionRecord
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_gks_record._CivicGksAssertionRecord
   :members:
   :show-inheritance:

CivicGksPredictiveAssertion
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_gks_record.CivicGksPredictiveAssertion
   :members:
   :show-inheritance:

CivicGksDiagnosticAssertion
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_gks_record.CivicGksDiagnosticAssertion
   :members:
   :show-inheritance:

CivicGksPrognosticAssertion
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_gks_record.CivicGksPrognosticAssertion
   :members:
   :show-inheritance:

CivicGksWriter
~~~~~~~~~~~~~~

.. autoclass:: civicpy.exports.civic_gks_writer.CivicGksWriter
   :members:

Examples
~~~~~~~~

Here's an example of how to export all assertions to GKS JSON::

    from civicpy import civic
    from civicpy.exports.civic_gks_writer import CivicGksWriter
    from civicpy.exports.civic_vcf_record import CivicVcfRecord

    records = []

    for assertion in civic.get_all_assertions():
      if assertion.is_valid_for_gks_json():
        try:
          if assertion.assertion_type == "DIAGNOSTIC":
            gks_record = CivicGksDiagnosticAssertion(assertion)
          elif assertion.assertion_type == "PREDICTIVE":
            gks_record = CivicGksPredictiveAssertion(assertion)
          else:
            gks_record = CivicGksPrognosticAssertion(assertion)
        except CivicGksRecordError:
          continue

        records.append(gks_record)
    CivicGksWriter("gks.json", records)

Here's an example of how to export all assertions approved by a specific organization that are
ready for submission to ClinVar.::

    from civicpy import civic
    from civicpy.exports.civic_gks_writer import CivicGksWriter
    from civicpy.exports.civic_vcf_record import CivicVcfRecord

    records = []
    organization_id = 1

    for assertion in civic.get_all_assertions_ready_for_clinvar_submission_for_org(organization_id):
      if assertion.is_valid_for_gks_json():
        try:
          if assertion.assertion_type == "DIAGNOSTIC":
            gks_record = CivicGksDiagnosticAssertion(assertion)
          elif assertion.assertion_type == "PREDICTIVE":
            gks_record = CivicGksPredictiveAssertion(assertion)
          else:
            gks_record = CivicGksPrognosticAssertion(assertion)
        except CivicGksRecordError:
          continue

        records.append(gks_record)
    CivicGksWriter("gks.json", records)
