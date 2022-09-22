default:
	@echo "Please specify an action: test | update_test_cache"
lint:
	@flake8 civicpy
test:
	@python -m pytest
update_test_cache:
	@python -m civicpy.cli update --hard --cache-save-path civicpy/data/test_cache.pkl