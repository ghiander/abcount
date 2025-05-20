import json
import logging

import rdkit
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
    def __init__(
        self,
        acid_defs_filepath=fps.acid_defs_filepath,
        base_defs_filepath=fps.base_defs_filepath,
    ):
        """
        Counts acidic and basic functional groups in a molecule using SMARTS pattern matching.

        Args:
            acid_defs_filepath (str): Path to custom acidic SMARTS definitions (optional).
            base_defs_filepath (str): Path to custom basic SMARTS definitions (optional).

        Attributes:
            acid_matcher (SmartsMatcherJson): Matcher for acidic functional group patterns.
            base_matcher (SmartsMatcherJson): Matcher for basic functional group patterns.
        """
        self.acid_defs_filepath = acid_defs_filepath
        self.base_defs_filepath = base_defs_filepath
        logger.debug(f"Loading Acidic SMARTS from path: {self.acid_defs_filepath}")
        logger.debug(f"Loading Basic SMARTS from path: {self.base_defs_filepath}")
        self.acid_matcher = SmartsMatcherJson(self.acid_defs_filepath)
        self.base_matcher = SmartsMatcherJson(self.base_defs_filepath)

    def count_acid_and_bases(self, mol: rdkit.Chem.Mol):
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
