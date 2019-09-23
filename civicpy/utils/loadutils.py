from civicpy import civic, TEST_CACHE_PATH, LOCAL_CACHE_PATH


def hard_update(local_cache_path=LOCAL_CACHE_PATH):
    """
    Construct a cache from the CIViC API and save to local_cache_path.

    :param local_cache_path:    A filepath destination string for storing the cache.
                                This parameter defaults to LOCAL_CACHE_PATH.

    """
    civic.update_cache(from_remote_cache=False, local_cache_path=local_cache_path)


if __name__ == '__main__':
    # This is for updating the test cache when executing this file as a script.
    hard_update(TEST_CACHE_PATH)
