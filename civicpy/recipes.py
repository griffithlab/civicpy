from civicpy import civic, exports
import requests


def get_variant_phenotypes_with_suggested_changes(variant_id_list):
    '''for each variant, yields evidence and associated phenotypes, both current and suggested'''
    variants = civic.get_variants_by_ids(variant_id_list)
    evidence = list()
    for variant in variants:
        evidence.extend(variant.evidence)
    for e in evidence:
        suggested_changes_url = f'https://civicdb.org/api/evidence_items/{e.id}/suggested_changes'
        resp = requests.get(suggested_changes_url)
        resp.raise_for_status()
        suggested_changes = dict()
        for suggested_change in resp.json():
            pheno_changes = suggested_change['suggested_changes'].get('phenotype_ids', None)
            if pheno_changes is None:
                continue
            a, b = pheno_changes
            added = set(b) - set(a)
            deleted = set(a) - set(b)
            rid = suggested_change['id']
            suggested_changes[rid] = {'added': added, 'deleted': deleted}
        yield e, {'suggested_changes': suggested_changes, 'current': set([x.id for x in e.phenotypes])}


def get_variant_phenotypes_with_suggested_changes_merged(variant_id_list):
    '''for each variant, yields evidence and merged phenotype from applying suggested changes to current'''
    for evidence, phenotype_status in get_variant_phenotypes_with_suggested_changes(variant_id_list):
        final = phenotype_status['current']
        for rid in sorted(phenotype_status['suggested_changes']):
            changes = phenotype_status['suggested_changes'][rid]
            final = final - changes['deleted']
            final = final | changes['added']
        if final:
            yield evidence, final


def export_all_variants_to_vcf(file):
    w = exports.VCFWriter(file)
    all_variants = civic.get_all_variants()
    w.addrecords(all_variants)
    w.writerecords()
