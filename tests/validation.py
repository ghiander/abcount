import csv
import os

import _match
import _report
from _match import ABValidator
from _match import MatchCounter
from _report import ReportEntryFactory
from _report import ReportGenerator
from rdkit import Chem
from rdkit.Chem import Mol

from abcount import ABCounter


def run_validation(
    mol: Mol,
    abcounter: ABCounter,
    validator: ABValidator,
    acid_counter: MatchCounter,
    base_counter: MatchCounter,
    acid_reporter: ReportGenerator,
    base_reporter: ReportGenerator,
):
    counts = abcounter.count_acid_and_bases(mol)

    # Acid evaluation
    predicted_acid = counts["acid"]
    expected_acid = int(row["pka_acid_num"])
    outcome = validator.generate_outcome(expected_acid, predicted_acid)
    acid_counter.count(outcome)

    # Acid FP/FN reporting
    entry = ReportEntryFactory.generate_entry(
        abcounter.acid_matcher, outcome, row["canonical_smiles"]
    )
    acid_reporter.add(entry)

    # Base evaluation
    predicted_base = counts["base"]
    expected_base = int(row["pka_base_num"])
    outcome = validator.generate_outcome(expected_base, predicted_base)
    base_counter.count(outcome)

    # Base FP/FN reporting
    entry = ReportEntryFactory.generate_entry(
        abcounter.base_matcher, outcome, row["canonical_smiles"]
    )
    base_reporter.add(entry)


# Load input data
file_path = os.path.join(os.path.dirname(__file__), "chembl_example_1.csv")
with open(file_path, mode="r", newline="") as file:
    reader = csv.DictReader(file)
    data = [row for row in reader]

# Initialise objects
abcounter = ABCounter()
validator = _match.ABValidator
a_counter = _match.MatchCounter()
b_counter = _match.MatchCounter()
a_reporter = _report.ReportGenerator()
b_reporter = _report.ReportGenerator()

# Evaluate predictions
for row in data:
    mol = row["mol"] = Chem.MolFromSmiles(row["canonical_smiles"])
    run_validation(
        mol, abcounter, validator, a_counter, b_counter, a_reporter, b_reporter
    )

print(f"Acidic groups — {a_counter.make_report()}")
print(f"Basic groups — {b_counter.make_report()}")


# Generate report
a_reporter.save_fp_report(prefix="acid_")
b_reporter.save_fp_report(prefix="base_")
a_reporter.save_fn_report(prefix="acid_")
b_reporter.save_fn_report(prefix="base_")
