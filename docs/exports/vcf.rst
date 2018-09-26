VCF
===

VCFs are written using the :class:`VCFWriter` class, to which you add :class:`Assertion`,
:class:`Variant`, or :class:`Evidence` records. The VCF output has one line per evidence.
Passing :class:`Assertion` or :class:`Variant` objects will expand the record to all
evidence describing it.

VCFWriter
---------

.. class:: VCFWriter
