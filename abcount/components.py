import json
import logging

from rdkit import Chem

from abcount.config import fps
from abcount.model import GroupTypeAttribute


logger = logging.getLogger(__name__)


class SmartsMatcher:
    """Plain SMARTS matcher with count functionality."""

    def __init__(self, definitions_list: list):
        self._load_smarts(definitions_list)

    def _load_smarts(self, definitions_list):
        # Loaded as a tuple to retain original SMARTS for debugging
        self.definitions_tup = [(Chem.MolFromSmarts(d), d) for d in definitions_list]

    def count_matches(self, mol):
        count = 0
        for t in self.definitions_tup:
            matches = len(mol.GetSubstructMatches(t[0]))
            if matches:
                logger.debug(f"SMARTS: {t[1]} - Matches: {matches}")
            count += matches
        return count


class SmartsMatcherJson(SmartsMatcher):
    """Adapter for SmartsMatcher."""

    def __init__(self, definitions_fp):
        self.definitions_fp = definitions_fp
        self._load_smarts()

    def _load_smarts(self):
        definitions_dict = self._load_json()
        # Loaded as a tuple to retain original SMARTS for debugging
        self.definitions_tup = [
            (Chem.MolFromSmarts(d["smarts"]), d) for d in definitions_dict
        ]

    def _load_json(self):
        with open(self.definitions_fp) as f:
            return json.load(f)


class ABCounter:
    def __init__(self):
        """
        Counts acidic and basic functional groups in a molecule using SMARTS pattern matching.

        Attributes:
            acid_matcher (SmartsMatcherJson): Matcher for acidic functional group patterns.
            base_matcher (SmartsMatcherJson): Matcher for basic functional group patterns.
        """

        self.acid_matcher = SmartsMatcherJson(fps.acid_defs_filepath)
        self.base_matcher = SmartsMatcherJson(fps.base_defs_filepath)

    def count_acid_and_bases(self, mol):
        """
        Count acidic and basic groups in the given molecule.

        Args:
            mol (rdkit.Chem.Mol): An RDKit molecule object.

        Returns:
            dict: A dictionary with counts keyed by GroupTypeAttribute.ACID and
                  GroupTypeAttribute.BASE.
        """
        acid_count = self.acid_matcher.count_matches(mol)
        base_count = self.base_matcher.count_matches(mol)
        return {
            GroupTypeAttribute.ACID: acid_count,
            GroupTypeAttribute.BASE: base_count,
        }


if __name__ == "__main__":
    # Test
    smiles = "OC(CC)(C)C#C"
    acid_exp = 0
    base_exp = 3
    mol = Chem.MolFromSmiles(smiles)
    abc = ABCounter()
    print(abc.count_acid_and_bases(mol))
