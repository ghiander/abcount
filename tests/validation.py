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
tp_acid = fp_acid = tn_acid = fn_acid = 0
tp_base = fp_base = tn_base = fn_base = 0
for row in data:
    counts = abcounter.count_acid_and_bases(row["mol"])

    # Acid evaluation
    predicted_acid = counts["acid"]
    expected_acid = int(row["pka_acid_num"])

    # Strict validation for 0 or 1 groups
    for n in (0, 1):
        if expected_acid == n:
            if n == 0 and predicted_acid == n:  # true negative
                tn_acid += 1
            if n == 1 and predicted_acid == n:  # true positive
                tp_acid += 1
            if predicted_acid > n:  # false positive
                fp_acid += 1
                print(row)
            if predicted_acid < n:  # only applies to n == 1
                fn_acid += 1
    # Relaxed validation in case of 2 groups
    if expected_acid == 2:
        if (
            predicted_acid >= expected_acid - 1
        ):  # true positive if at least 1 group is matched
            tp_acid += 1
        else:
            fn_acid += 1

    # Base evaluation
    predicted_base = counts["base"]
    expected_base = int(row["pka_base_num"])

    # Strict validation for 0 or 1 groups
    for n in (0, 1):
        if expected_base == n:
            if n == 0 and predicted_base == n:  # true negative
                tn_base += 1
            if n == 1 and predicted_base == n:  # true positive
                tp_base += 1
            if predicted_base > n:  # false positive
                fp_base += 1
            if predicted_base < n:  # only applies to n == 1
                fn_base += 1
    # Relaxed validation in case of 2 groups
    if expected_base == 2:
        if (
            predicted_base >= expected_acid - 1
        ):  # true positive if at least 1 group is matched
            tp_base += 1
        else:
            fn_base += 1

print(f"Acidic groups — TP: {tp_acid} | FP: {fp_acid} | TN: {tn_acid} | FN: {fn_acid}")
print(f"Basic groups — TP: {tp_base} | FP: {fp_base} | TN: {tn_base} | FN: {fn_base}")
