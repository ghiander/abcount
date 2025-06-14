.PHONY: build install upload upload_test test coverage create_definitions

definitions:
	python abcount/_definitions.py

build:
	uv build

test:
	pytest -vss tests/test.py
	pytest -vss tests/test_validation.py

validate:
	cd tests && python validation.py

install:
	pip install dist/abcount-*.tar.gz --force-reinstall

upload_test:
	uv publish --index testpypi

upload:
	twine upload -r pypi dist/*

coverage:
	pytest --cov=src test/ -vv --cov-report term-missing
