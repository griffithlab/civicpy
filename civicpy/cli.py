import click
from civicpy import LOCAL_CACHE_PATH, civic
from civicpy.exports import VCFWriter


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('--soft/--hard', default=True,
              help='Hard-update from live API (slow) or \
              soft-update from daily precache (fast; default)')
@click.option('--cache-save-path',
              help=f'Filepath to save cache to. Default: {LOCAL_CACHE_PATH}',
              default=LOCAL_CACHE_PATH)
def update(soft, cache_save_path):
    """Updates CIViC content from server and stores to local cache file"""
    civic.update_cache(from_remote_cache=soft, local_cache_path=cache_save_path)

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-v', '--vcf-file-path', required=True,
              help="The file path to write the VCF to.")
@click.option('-i', '--include-status', required=True, multiple=True, type=click.Choice(['accepted', 'submitted', 'rejected']),
              help="Limits the variants and annotations in the VCF to only the ones that match the given statuses. \
              May be specified more than once.")
def create_vcf(vcf_file_path, include_status):
    """Create a VCF file of CIViC variants"""
    with open(vcf_file_path, "w") as fh:
        writer = VCFWriter(fh)
        for variant in civic.get_all_variants(include_status=include_status):
            if variant.is_valid_for_vcf():
                writer.addrecord(variant)
        writer.writerecords()

if __name__ == '__main__':
    cli()
