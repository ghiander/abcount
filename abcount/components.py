import json

from rdkit import Chem

from abcount.config import fps
from abcount.model import GroupTypeAttribute


# class SmartsMatcher:
#     def __init__(self, definitions_fp):
#         self.definitions_fp = definitions_fp
#         self.definitions = self._load_smarts_obj()

#     def _load_smarts_obj(self):
#         defs_dict = self._load_json()
#         return [Chem.MolFromSmarts(d["smarts"]) for d in defs_dict]

#     def _load_json(self):
#         with open(self.definitions_fp) as f:
#             return json.load(f)

#     def count_matches(self, mol):
#         count = 0
#         for s in self.definitions:
#             count += len(mol.GetSubstructMatches(s))
#         return count


class SmartsMatcher:
    """Plain SMARTS matcher with count functionality."""

    def __init__(self, definitions_list: list):
        self.definitions = self._load_smarts_obj(definitions_list)

    def _load_smarts_obj(self, definitions_list):
        return [Chem.MolFromSmarts(d) for d in definitions_list]

    def count_matches(self, mol):
        count = 0
        for s in self.definitions:
            count += len(mol.GetSubstructMatches(s))
        return count


class SmartsMatcherJson(SmartsMatcher):
    """Adapter for SmartsMatcher."""

    def __init__(self, definitions_fp):
        self.definitions_fp = definitions_fp
        self.definitions = self._load_smarts_obj()

    def _load_smarts_obj(self):
        defs_dict = self._load_json()
        return [Chem.MolFromSmarts(d["smarts"]) for d in defs_dict]

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
