import click
from civicpy import LOCAL_CACHE_PATH, civic


@click.group()
def cli():
    pass


@cli.command()
@click.option('--soft/--hard', default=True,
              help='Hard-update from live API (slow) or \
              soft-update from daily precache (fast; default)')
@click.option('--cache-save-path',
              help=f'Filepath to save cache to. Default: {LOCAL_CACHE_PATH}',
              default=LOCAL_CACHE_PATH)
def update(soft, cache_save_path):
    """Updates CIViC content from server and stores to local cache file"""
    civic.update_cache(from_remote_cache=soft, local_cache_path=cache_save_path)


if __name__ == '__main__':
    cli()
