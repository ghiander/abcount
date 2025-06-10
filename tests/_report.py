from dataclasses import dataclass
from datetime import datetime

import pandas as pd
from _match import FalseNegative
from _match import FalsePositive
from _match import PredictionOutcome
from rdkit import Chem

from abcount.components import SmartsMatcher


class ReportEntry:
    pass


class NullReportEntry(ReportEntry):
    pass


@dataclass
class FPReportEntry(ReportEntry):
    target_smiles: str
    matches_list: list
    expected_matches: int


@dataclass
class FNReportEntry(ReportEntry):
    target_smiles: str
    expected_matches: int


class ReportEntryFactory:
    @staticmethod
    def generate_entry(
        matcher: SmartsMatcher, outcome: PredictionOutcome, target_smiles: str
    ):
        if isinstance(outcome, FalsePositive):
            mol = Chem.MolFromSmiles(target_smiles)
            if not mol:
                raise RuntimeError(
                    "`outcome` was provided but `target_smiles` are not valid."
                )
            matches_list = matcher.generate_matches_list(mol)
            return FPReportEntry(
                target_smiles=target_smiles,
                matches_list=matches_list,
                expected_matches=outcome.expected_matches,
            )
        if isinstance(outcome, FalseNegative):
            return FNReportEntry(
                target_smiles=target_smiles, expected_matches=outcome.expected_matches
            )
        else:
            # Just to avoid returning Nones at lower level as good old Bob recommends
            return NullReportEntry()


class ReportGenerator:
    def __init__(self):
        self._fp_columns = ["target", "smarts", "matches", "total_expected_matches"]
        self.fp_df = pd.DataFrame(columns=self._fp_columns)
        self._fn_columns = ["target", "total_expected_matches"]
        self.fn_df = pd.DataFrame(columns=self._fn_columns)

    def add(self, entry: ReportEntry):
        if isinstance(entry, FPReportEntry):
            self._add_fp(entry)
        elif isinstance(entry, FNReportEntry):
            self._add_fn(entry)
        elif not isinstance(entry, ReportEntry):
            raise TypeError("`ReportEntry` instance is required.")

    def _add_fp(self, entry: FPReportEntry):
        target_smiles = entry.target_smiles
        expected_matches = entry.expected_matches
        for e in entry.matches_list:
            ungrouped_entry = (target_smiles, e[0]["smarts"], e[1], expected_matches)
            self.fp_df.loc[len(self.fp_df)] = ungrouped_entry

    def _add_fn(self, entry: FNReportEntry):
        target_smiles = entry.target_smiles
        expected_matches = entry.expected_matches
        self.fn_df.loc[len(self.fn_df)] = [target_smiles, expected_matches]

    def save_fp_report(self, prefix: str = ""):
        self.fp_df.to_csv(
            f"{prefix}fp_report_{ReportGenerator._get_timestamp()}.csv", index=False
        )

    def save_fn_report(self, prefix: str = ""):
        self.fn_df.to_csv(
            f"{prefix}fn_report_{ReportGenerator._get_timestamp()}.csv", index=False
        )

    @staticmethod
    def _get_timestamp():
        return datetime.utcnow().strftime("%Y%m%dT%H%M%S")
