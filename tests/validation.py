import csv
import os

from rdkit import Chem

from abcount import ABCounter

# Load input data
file_path = os.path.join(os.path.dirname(__file__), "chembl_example_1.csv")
with open(file_path, mode="r", newline="") as file:
    reader = csv.DictReader(file)
    data = [row for row in reader]
for row in data:
    row["mol"] = Chem.MolFromSmiles(row["canonical_smiles"])


# Evaluate predictions
abcounter = ABCounter()
correct_acid = incorrect_acid = 0
correct_base = incorrect_base = 0
for row in data:
    counts = abcounter.count_acid_and_bases(row["mol"])

    # Acid evaluation
    predicted_acid = counts["acid"]
    expected_acid = int(row["pka_acid_num"])

    # Strict validation in case of no groups
    if expected_acid == 0:
        if predicted_acid == 0:
            correct_acid += 1
        else:
            incorrect_acid += 1
    # Relaxed validation in case of groups
    elif expected_acid in (1, 2):
        if predicted_acid >= expected_acid:
            correct_acid += 1
        else:
            incorrect_acid += 1

    # Base evaluation
    predicted_base = counts["base"]
    expected_base = int(row["pka_base_num"])

    # Strict validation in case of no groups
    if expected_base == 0:
        if predicted_base == 0:
            correct_base += 1
        else:
            incorrect_base += 1
    # Relaxed validation in case of groups
    elif expected_base in (1, 2):
        if predicted_base >= expected_base:
            correct_base += 1
        else:
            incorrect_base += 1


print("Acidic groups — Correct:", correct_acid, "| Incorrect:", incorrect_acid)
print("Basic groups — Correct:", correct_base, "| Incorrect:", incorrect_base)
