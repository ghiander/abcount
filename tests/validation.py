import csv
import os

import _match
import _report
from rdkit import Chem

from abcount import ABCounter

# Load input data
file_path = os.path.join(os.path.dirname(__file__), "chembl_example_1.csv")
with open(file_path, mode="r", newline="") as file:
    reader = csv.DictReader(file)
    data = [row for row in reader]
for row in data:
    row["mol"] = Chem.MolFromSmiles(row["canonical_smiles"])


validator = _match.ABValidator
a_counter = _match.MatchCounter()
b_counter = _match.MatchCounter()

entry_factory = _report.ReportEntryFactory
a_reporter = _report.ReportGenerator()
b_reporter = _report.ReportGenerator()

# Evaluate predictions
abcounter = ABCounter()
tp_acid = fp_acid = tn_acid = fn_acid = 0
tp_base = fp_base = tn_base = fn_base = 0
for row in data:
    counts = abcounter.count_acid_and_bases(row["mol"])

    # Acid evaluation
    predicted_acid = counts["acid"]
    expected_acid = int(row["pka_acid_num"])
    outcome = validator.generate_outcome(expected_acid, predicted_acid)
    a_counter.count(outcome)

    # Acid FP/FN reporting
    entry = entry_factory.generate_entry(
        abcounter.acid_matcher, outcome, row["canonical_smiles"]
    )
    a_reporter.add(entry)

    # Base evaluation
    predicted_base = counts["base"]
    expected_base = int(row["pka_base_num"])
    outcome = validator.generate_outcome(expected_base, predicted_base)
    b_counter.count(outcome)

    # Base FP/FN reporting
    entry = entry_factory.generate_entry(
        abcounter.base_matcher, outcome, row["canonical_smiles"]
    )
    b_reporter.add(entry)


print(f"Acidic groups — {a_counter.make_report()}")
print(f"Basic groups — {b_counter.make_report()}")


# Generate report
a_reporter.save_fp_report(prefix="acid_")
b_reporter.save_fp_report(prefix="base_")
a_reporter.save_fn_report(prefix="acid_")
b_reporter.save_fn_report(prefix="base_")
