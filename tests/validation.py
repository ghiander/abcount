import csv
import os

import _match
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

    # Base evaluation
    predicted_base = counts["base"]
    expected_base = int(row["pka_base_num"])
    outcome = validator.generate_outcome(expected_base, predicted_base)
    b_counter.count(outcome)

print(f"Acidic groups — {a_counter.make_report()}")
print(f"Basic groups — {b_counter.make_report()}")
