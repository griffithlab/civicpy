.. module:: civicpy.exports

Exporting
=========

CIViCpy supports exporting of CIViC records to Variant Call Format (VCF) files.
This enables downstream analyses such as integrating with `IGV`_, `VEP`_, and
other common bioinformatics tools. VCF exports are maintained via the :mod:`exports`
module.

Other file formats are planned for future releases. Suggestions are welcome on our
`GitHub issues page <https://github.com/griffithlab/civicpy/issues>`_.

VCF
---

VCFs are written using the :class:`VCFWriter` class, to which you add :class:`Assertion`,
:class:`Variant`, or :class:`Evidence` records. The VCF output has one line per evidence.
Passing :class:`Assertion` or :class:`Variant` objects will expand the record to all
evidence describing it.

VCFWriter
~~~~~~~~~

.. class:: VCFWriter

Example
~~~~~~~~~

Here's an example of how to export all variants from CIViC to VCF::

	from civicpy import civic, exports

	with open('civic_variants.vcf', 'w', newline='') as file:
		w = exports.VCFWriter(file)
		all_variants = civic.get_all_variants()
		w.addrecords(all_variants)
		w.writerecords()

.. _`IGV`: https://software.broadinstitute.org/software/igv/
.. _`VEP`: https://useast.ensembl.org/info/docs/tools/vep/index.html