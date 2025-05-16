.PHONY: build install upload upload_test test coverage

build:
	uv build

install:
	pip install dist/abcount-*.tar.gz --force-reinstall

upload_test:
	twine upload -r testpypi dist/*

upload:
	twine upload -r pypi dist/*

test:
	pytest -vss tests/test.py

coverage:
	pytest --cov=src test/ -vv --cov-report term-missing
