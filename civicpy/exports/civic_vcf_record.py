import vcfpy
import requests
import re
from civicpy.civic import GeneVariant

class CivicVcfRecord(vcfpy.Record):
    """
    :param civic.GeneVariant variant: A :class:`civic.GeneVariant` object to convert to a :class:`civic.exports.CivicVcfRecord` object (inherits from vcfpy.Record)
    """

    def __init__(self, variant, include_status=["submitted", "accepted"]):
        if not isinstance(variant, GeneVariant):
            raise Exception('Variant is not a GeneVariant.')
        self.variant = variant

        if not self.variant.is_valid_for_vcf(emit_warnings=True):
            raise Exception('Variant is not valid for VCF.')

        if self.vcf_coordinates() is None:
            raise Exception("Variant doesn't have any coordinates.")

        (start, ref, alt) = self.vcf_coordinates()
        info = vcfpy.OrderedDict({
            ('GN', self.variant.gene.name),
            ('VT', self.sanitized_name()),
            ('CSQ', tuple(self.csq(include_status))),
        })

        super().__init__(
            CHROM=self.variant.coordinates.chromosome,
            POS=start,
            ID=[str(self.variant.id)],
            REF=ref,
            ALT=[self.alt_record(alt)],
            QUAL=None,
            FILTER=['.'],
            INFO=info,
        )

    def vcf_coordinates(self):
        ensembl_server = "https://grch37.rest.ensembl.org"
        if self.variant.coordinates.reference_build != 'GRCh37':
            return
        if self.variant.is_insertion:
            if not self.variant.coordinates.representative_transcript:
                return
            else:
                start = self.variant.coordinates.start
                ext = "/sequence/region/human/{}:{}-{}".format(self.variant.coordinates.chromosome, start, start)
                r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "text/plain"})
                r.raise_for_status()
                if self.variant.coordinates.reference_bases == None or self.variant.coordinates.reference_bases == '-' or self.variant.coordinates.reference_bases == '':
                    ref = r.text
                else:
                    ref = "{}{}".format(r.text, self.variant.coordinates.reference_bases)
                alt = "{}{}".format(r.text, self.variant.coordinates.variant_bases)
        elif self.variant.is_deletion:
            if not self.variant.coordinates.representative_transcript:
                return
            else:
                start = self.variant.coordinates.start - 1
                ext = "/sequence/region/human/{}:{}-{}".format(self.variant.coordinates.chromosome, start, start)
                r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "text/plain"})
                r.raise_for_status()
                ref = "{}{}".format(r.text, self.variant.coordinates.reference_bases)
                if self.variant.coordinates.variant_bases == None or self.variant.coordinates.variant_bases == '-' or self.variant.coordinates.variant_bases == '':
                    alt = r.text
                else:
                    alt = "{}{}".format(r.text, self.variant.coordinates.variant_bases)
        else:
            start = self.variant.coordinates.start
            ref = self.variant.coordinates.reference_bases
            alt = self.variant.coordinates.variant_bases
        return (start, ref, alt)

    def alt_record(self, alt):
        if self.variant.is_insertion:
            return vcfpy.Substitution(type_="INS", value=alt)
        elif self.variant.is_deletion:
            return vcfpy.Substitution(type_="DEL", value=alt)
        else:
            if len(alt) == 1:
                return vcfpy.Substitution(type_="SNV", value=alt)
            else:
                return vcfpy.Substitution(type_="MNV", value=alt)

    def csq_alt(self):
        if self.variant.coordinates.reference_build != 'GRCh37':
            return
        if self.variant.is_insertion:
            if not self.variant.coordinates.representative_transcript:
                return
            else:
                return self.variant.coordinates.variant_bases
        elif self.variant.is_deletion:
            if not self.variant.coordinates.representative_transcript:
                return
            else:
                return "-"
        else:
            return self.variant.coordinates.variant_bases

    def hgvs_c(self):
        if self.variant.coordinates.representative_transcript:
            hgvs_cs = [e for e in self.variant.hgvs_expressions if (':c.' in e) and (self.variant.coordinates.representative_transcript in e)]
            return hgvs_cs[0] if len(hgvs_cs) == 1 else ''
        else:
            return ''

    def hgvs_p(self):
        if self.variant.coordinates.representative_transcript:
            hgvs_ps = [e for e in self.variant.hgvs_expressions if (':p.' in e) and (self.variant.coordinates.representative_transcript in e)]
            return hgvs_ps[0] if len(hgvs_ps) == 1 else ''
        else:
            return ''

    def sanitized_name(self):
        name = self.variant.name
        regex = re.compile(r"^([A-Z]+)([0-9]+)(=)(.*)$")
        match = regex.match(name)
        if match is not None:
            name = "".join([match.group(1), match.group(2), match.group(1), match.group(4)])
        return name

    def csq(self, include_status=None):
        if self.csq_alt() is None:
            return []
        else:
            csq = []
            for mp in self.variant.molecular_profiles:
                for evidence in mp.evidence:
                    if include_status is not None and evidence.status not in include_status:
                        continue
                    csq.append('|'.join([
                        self.csq_alt(),
                        '&'.join(map(lambda t: t.name, self.variant.variant_types)),
                        self.variant.gene.name,
                        str(self.variant.gene.entrez_id),
                        'transcript',
                        str(self.variant.coordinates.representative_transcript),
                        self.hgvs_c(),
                        self.hgvs_p(),
                        self.sanitized_name(),
                        str(self.variant.id),
                        '&'.join(self.variant.variant_aliases),
                        "https://civicdb.org/links/variants/{}".format(self.variant.id),
                        mp.sanitized_name(),
                        str(mp.id),
                        '&'.join(mp.aliases),
                        "https://civicdb.org/links/molecular-profiles/{}".format(mp.id),
                        '&'.join(self.variant.hgvs_expressions),
                        str(self.variant.allele_registry_id),
                        '&'.join(self.variant.clinvar_entries),
                        str(mp.molecular_profile_score),
                        "evidence",
                        str(evidence.id),
                        "https://civicdb.org/links/evidence/{}".format(evidence.id),
                        "{} ({})".format(evidence.source.citation_id, evidence.source.source_type),
                        str(evidence.variant_origin),
                        evidence.status,
                        str(evidence.significance or ''),
                        str(evidence.evidence_direction or ''),
                        evidence.disease.name if evidence.disease is not None else "",
                        '&'.join([str(therapy) for therapy in evidence.therapies]),
                        str(evidence.therapy_interaction_type or ""),
                        '&'.join(["{} (HPO ID {})".format(phenotype.name, phenotype.hpo_id) for phenotype in evidence.phenotypes]),
                        evidence.evidence_level,
                        str(evidence.rating),
                        "",
                        "",
                        "",
                        "",
                        "",
                    ]))
                for assertion in mp.assertions:
                    if include_status is not None and assertion.status not in include_status:
                        continue
                    csq.append('|'.join([
                        self.csq_alt(),
                        '&'.join(map(lambda t: t.name, self.variant.variant_types)),
                        self.variant.gene.name,
                        str(self.variant.gene.entrez_id),
                        'transcript',
                        str(self.variant.coordinates.representative_transcript),
                        self.hgvs_c(),
                        self.hgvs_p(),
                        self.sanitized_name(),
                        str(self.variant.id),
                        '&'.join(self.variant.variant_aliases),
                        "https://civicdb.org/links/variants/{}".format(self.variant.id),
                        mp.sanitized_name(),
                        str(mp.id),
                        '&'.join(mp.aliases),
                        "https://civicdb.org/links/molecular-profiles/{}".format(mp.id),
                        '&'.join(self.variant.hgvs_expressions),
                        str(self.variant.allele_registry_id),
                        '&'.join(self.variant.clinvar_entries),
                        str(mp.molecular_profile_score),
                        "assertion",
                        str(assertion.id),
                        "https://civicdb.org/links/assertion/{}".format(assertion.id),
                        "",
                        str(assertion.variant_origin),
                        assertion.status,
                        assertion.significance,
                        assertion.assertion_direction,
                        str(assertion.disease),
                        '&'.join([str(therapy) for therapy in assertion.therapies]),
                        str(assertion.therapy_interaction_type or ''),
                        "",
                        "",
                        "",
                        "&".join([acmg_code.code for acmg_code in assertion.acmg_codes]),
                        str(assertion.amp_level or ''),
                        assertion.format_nccn_guideline(),
                        str(assertion.fda_regulatory_approval or ''),
                        str(assertion.fda_companion_test or ''),
                    ]))
            return csq
