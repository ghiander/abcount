# ABCount - AcidBaseCount
ABCount is a SMARTS-based tool that determines the number of acidic and basic groups in molecules.

## How to install the tool
ABCount can be installed from pypi (https://pypi.org/project/abcount).
```bash
pip install abcount
```

## Usage
```python
from rdkit import Chem
from abcount import ABCounter

mol = Chem.MolFromSmiles("[nH]1nnnc1-c3c2[nH]ncc2ccc3")
abc = ABCounter()
abc.count_acid_and_bases(mol)
```
```python
{'acid': 2, 'base': 2}
```

## SMARTS definitions Source
The SMARTS patterns used in this project were obtained from the following sources. Note that definitions are not deduplicated, hence require curation to avoid redundant matching.

* Pan, X.; Wang, H.; Li, C.; Zhang, J. Z. H.; Ji, C., **MolGpka: A Web Server for Small Molecule pKa Prediction Using a Graph-Convolutional Neural Network**
*Journal of Chemical Information and Modeling* **2021**, *61* (7), 3159–3165.
DOI: [10.1021/acs.jcim.1c00075](https://doi.org/10.1021/acs.jcim.1c00075)
* Some manually curated definitions.

## Some useful commands
- Generate acidic and basic definitions from aggregated data: `python abcount/_definitions.py`. A follow up on how definitions can be curated will be provided.
- Run tests: `pytest -vss tests/test.py`
- Run validation: `python tests/validation.py`

## For developers
- The package was created using `uv` (https://docs.astral.sh/uv/).
- The package can be installed from the wheel in the `dist/` folder. When a new version needs to be released, a new wheel must be built. That can be done by changing the version of the package inside `setup.py` then calling `python setup.py bdist_wheel` and `python setup.py sdist` which will create a new build.
- The code can be automatically tested using `pytest -vss tests/test.py` which requires `pytest` to be installed.
- The `Makefile` can also be used for building (`make build`) or testing (`make test`).
- Before committing new code, please always check that the style and syntax are compliant using `pre-commit`.

### Setting up your development environment
The `pyproject.toml` already contains the optional dependencies needed for development. Follow these steps to set up the environment.
```bash
# Make sure you have got Python >= 3.10
python --version
> Python 3.10.16

# Installs `abcount` in editable mode and with dev dependencies
pip install -e .[dev]
> ...
> Successfully installed abcount ...

# Setup pre-commit hooks
pre-commit install
> pre-commit installed at .git/hooks/pre-commit
```
