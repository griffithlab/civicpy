Usage
=====

CIViCpy is primarily designed to enable exploration of the content of CIViC through Python :class:`CivicRecord` objects.
:class:`Gene`, :class:`Variant`, :class:`Assertion`, :class:`Source`, and :class:`Evidence`
are all subclasses of :class:`CivicRecord`. While these object can be generated locally,

The CivicRecord classes
-----------------------

.. class:: CivicRecord

	As a base class, :class:`CivicRecord` is used to define the characteristic of all records in CIViC. This class is not
	intended to be invoked directly by the end user, but provided for documentation of shared methods and variables in
	child classes.

	.. method:: update(allow_partial=True, force=False, **kwargs)

		Updates the record object from the cache or the server. The `allow_partial` flag will
		update the record according to the contents of CACHE, without requiring all attributes to be assigned. The
		`force` flag is used to force an update from the server, even if a full record exists in the cache. Keyword
		arguments may be passed to `kwargs`, which will update the corresponding attributes of the
		:class:`CivicRecord` instance.

	.. method:: site_link()

		Returns a URL to the record on the CIViC web application.

.. class:: Gene

.. class:: Variant

.. class:: Evidence

.. class:: Assertion

.. class:: Source