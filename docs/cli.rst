Commandline Utilities
=====================

Updating Your Local Cache
-------------------------

.. program-output:: civicpy update -h

Creating a VCF of CIViC Variants
--------------------------------

.. program-output:: civicpy create-vcf -h

Annotating a VCF with data from CIViC Variants
----------------------------------------------

.. program-output:: civicpy annotate-vcf -h

Create dereferenced GKS JSON for ClinVar
----------------------------------------

GKS export uses the Variation Normalizer REST service at
``http://127.0.0.1:8000/variation`` by default. Set the endpoint before running the
command when a different service should be used::

    export CIVICPY_VARIATION_NORMALIZER_URL=https://variation-normalizer.example/variation

.. program-output:: civicpy create-gks-json -h

Create referenced GKS Bundle JSON
---------------------------------

.. program-output:: civicpy create-gks-bundle -h
