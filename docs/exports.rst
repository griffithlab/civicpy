.. module:: exports

The **exports** module
======================

CIViCpy supports exporting of CIViC records to Variant Call Format (VCF) files.
This enables downstream analyses such as integrating with `IGV`_, `VEP`_, and
other common bioinformatics tools. VCF exports are maintained via the :mod:`exports`
module::

	>>>from civicpy import exports

Other file formats are planned for future releases. Suggestions are welcome on our
`GitHub issues page <https://github.com/griffithlab/civicpy/issues>`_.

VCF
---

VCFs are written using the :class:`VCFWriter` class, to which you add :class:`civic.Assertion`,
:class:`civic.Variant`, :class:`civic.Gene`, or :class:`civic.Evidence` records using the
``addrecord`` or ``addrecords`` methods. The VCF output
has one line per variant. Passing :class:`civic.Assertion`, :class:`civic.Gene`, or
:class:`civic.Evidence` objects will expand the record to all variants linked to those objects.

The ``addrecord`` function supports various variant types, depending on the curated coordinates
available for a variant. All variants require the chromosome name and start position. For SNVs
and complex variants the reference sequence and variant sequence information also need to be
available. By contrast, insertions require only variant sequence information and deletions
require only reference sequence information. Variants that do not meet these minimum requirements
will not be added and a warning message is emitted instead. Fusions and other variants with a
second set of coordinates are currently not supported.

In order to verify whether a variant can be added to a VCFWriter object, the convenience method
``is_valid_for_vcf`` can be called on a :class:`civic.Variant` object before calling
``addrecord``. Those variants that are unable to be exported into the VCF format are still
retrievable as CIViCpy records. Once all desired variants are added to the VCFWriter object,
``writerecords`` needs to be called to write the VCF file.

The variants added to the VCFWriter object are written to the VCF file, one VCF record for each
:class:`civic.Variant` object. If two variants share the same chromosome, start position, and
reference allele(s), they will not be combined into one VCF record but will instead be written
as separate VCF records. Additional CIViC data are added to the VCF as annotations to the
``CSQ`` (consequence) ``INFO`` field. All CIViC molcular profiles that this variant is a part of
are retrieved and all evidence items and assertions linked to these molecular profiles are added
to the CSQ field with one CSQ entry for each evidence item and/or assertion. Whether
a specific CSQ entry reflects an evidence item or an assertion is determined by the
``CIViC Entity Type`` CSQ field. To differentiate special characters in the field values from
field delimiters, spaces are replaced with underscores and other special characters are
hex-encoded. By utilizing the CSQ field for annotations, the resulting VCF is compatible for
import into Google BigQuery (git.io/bigquery-variant-annotation).

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

VCFWriter
~~~~~~~~~

.. autoclass:: exports.VCFWriter
   :members:

Example
~~~~~~~

Here's an example of how to export all variants from CIViC to VCF::

	from civicpy import civic, exports

	with open('civic_variants.vcf', 'w', newline='') as file:
		w = exports.VCFWriter(file)
		all_variants = civic.get_all_variants()
		w.addrecords(all_variants)
		w.writerecords()

.. _`IGV`: https://software.broadinstitute.org/software/igv/
.. _`VEP`: https://useast.ensembl.org/info/docs/tools/vep/index.html
