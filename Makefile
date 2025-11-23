.PHONY: build install upload upload_test test coverage create_definitions

definitions:
	python abcount/_definitions.py

build:
	uv build

test:
	pytest -vss tests/test_counter.py
	pytest -vss tests/test_classifier.py
	pytest -vss tests/test_ions.py
	pytest -vss tests/test_validation.py

create_git_tag:
	version=$$(grep -Po "version\s*=\s*['\"]\K[^'\"]+" pyproject.toml); \
	git tag v$$version
	git push --tags

validate:
	cd tests && python validation.py

install:
	pip install dist/abcount-*.tar.gz --force-reinstall

upload_test:
	version=$$(grep -Po "version\s*=\s*['\"]\K[^'\"]+" pyproject.toml); \
	uvx uv-publish --repository testpypi dist/abcount-$$version*

upload:
	version=$$(grep -Po "version\s*=\s*['\"]\K[^'\"]+" pyproject.toml); \
	twine upload -r pypi dist/abcount-$$version*

coverage:
	pytest --cov=src test/ -vv --cov-report term-missing
