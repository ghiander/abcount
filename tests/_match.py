class PredictionOutcome:
    pass


class TruePositive(PredictionOutcome):
    pass


class TrueNegative(PredictionOutcome):
    pass


class FalsePositive(PredictionOutcome):
    pass


class FalseNegative(PredictionOutcome):
    pass


class ABValidator:
    """Bespoke validator for AB classification."""

    @staticmethod
    def generate_outcome(expected_groups: int, predicted_groups: int):
        # Relaxed validation in case of 2 groups
        if expected_groups == 2:
            if (
                predicted_groups >= expected_groups - 1
            ):  # true positive if at least 1 group is matched
                return TruePositive()
            else:
                return FalseNegative()

        # Strict validation for 0 or 1 groups
        elif expected_groups in (0, 1):
            if (
                expected_groups == 0 and predicted_groups == expected_groups
            ):  # true negative
                return TrueNegative()
            if (
                expected_groups == 1 and predicted_groups == expected_groups
            ):  # true positive
                return TruePositive()
            if predicted_groups > expected_groups:  # false positive
                return FalsePositive()
            if predicted_groups < expected_groups:  # only applies to n == 1
                return FalseNegative()


class MatchCounter:
    def __init__(self):
        self.tps = 0
        self.tns = 0
        self.fps = 0
        self.fns = 0

    def count(self, match: PredictionOutcome):
        if isinstance(match, TruePositive):
            self.tps += 1
        elif isinstance(match, TrueNegative):
            self.tns += 1
        elif isinstance(match, FalsePositive):
            self.fps += 1
        elif isinstance(match, FalseNegative):
            self.fns += 1

    def make_report(self):
        return f"TP: {self.tps} | FP: {self.fps} | TN: {self.tns} | FN: {self.fns}"
