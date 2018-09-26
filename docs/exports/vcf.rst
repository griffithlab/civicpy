VCF
===

VCFs are written using the :class:`VCFWriter` class, to which you add :class:`Assertion`,
:class:`Variant`, or :class:`Evidence` records. The VCF output has one line per evidence.
Passing :class:`Assertion` or :class:`Variant` objects will expand the record to all
evidence describing it.

VCFWriter
---------

.. class:: VCFWriter

Example
-------

Here's an example of how to export all variants from CIViC to VCF::

	from civicpy import civic, exports

	with open('civic_variants.vcf', 'w', newline='') as file:
		w = exports.VCFWriter(file)
		all_variants = civic.get_all_variants()
		w.addrecords(all_variants)
		w.writerecords()
