import logging

from abcount.model.common import PKaAttribute
from abcount.model.isoelectric import NetCharge
from abcount.model.isoelectric import pH

logger = logging.getLogger(__name__)


class pIPredictor:
    @staticmethod
    def predict_input(
        pka_dict: dict, pka_attribute_cls=PKaAttribute, rounding_digits=2
    ):
        result = pIPredictor._predict_exact(
            pka_dict, pka_attribute_cls=pka_attribute_cls
        )
        return round(result, rounding_digits)

    @staticmethod
    def _predict_exact(pka_dict, pka_attribute_cls=PKaAttribute):
        # TODO: Bisect only works when charges cross zero.
        # Otherwise the upper or lower bounds may be used?
        # Investigate...
        return pIPredictor._predict_bisect(
            pka_dict, pka_attribute_cls=pka_attribute_cls
        )

    @staticmethod
    def _predict_bisect(pka_dict: dict, pka_attribute_cls=PKaAttribute):
        """Compute pI by finding pH where net charge = 0 using bisection.
        https://en.wikipedia.org/wiki/Bisection_method"""
        low, high = 0.0, 14.0

        q_low = NetChargeCalculator(pka_attribute_cls).calculate_at_pH(
            pka_dict, pH=pH(low)
        )
        q_high = NetChargeCalculator(pka_attribute_cls).calculate_at_pH(
            pka_dict, pH=pH(high)
        )
        if q_low * q_high > 0:
            # minus by plus must be minus
            raise ValueError(f"Net charge does not cross zero in [{low}, {high}]")

        while (high - low) > 1e-4:  # tolerance for iteration
            mid = (low + high) / 2
            q_mid = NetChargeCalculator(pka_attribute_cls).calculate_at_pH(
                pka_dict, pH=pH(mid)
            )
            logger.debug(
                {"low pH": low, "mid pH": mid, "high pH": high, "q_mid": q_mid}
            )

            if q_mid == 0:
                return mid

            # apply the bisect
            # negative charge - reduce higher bound
            elif q_low * q_mid < 0:
                high = mid
                q_high = q_mid
            # positive charge - reduce lower bound
            else:
                low = mid
                q_low = q_mid
        return 0.5 * (low + high)


class NetChargeCalculator:
    def __init__(self, pka_attribute_cls=PKaAttribute):
        self.__pka_attribute_cls__ = pka_attribute_cls
        self._set_attribute_lists()

    def _set_attribute_lists(self) -> None:
        self.acid_attribute_list = [
            self.__pka_attribute_cls__.ACID_1,
            self.__pka_attribute_cls__.ACID_2,
        ]
        self.base_attribute_list = [
            self.__pka_attribute_cls__.BASE_1,
            self.__pka_attribute_cls__.BASE_2,
        ]

    def calculate_at_pH(self, pka_dict: dict, pH: pH):
        nq = NetCharge()

        # Acidic charges
        for att in self.acid_attribute_list:
            pka_value = pka_dict[att]
            if pka_value:
                nq.add(ChargeCalculator.calculate_acid_charge(pka_value, pH.value))

        # Basic charges
        for att in self.base_attribute_list:
            pka_value = pka_dict[att]
            if pka_value:
                nq.add(ChargeCalculator.calculate_base_charge(pka_value, pH.value))
        return nq.value


class ChargeCalculator:
    """Calculates the fraction of protonated/deprotonated species based on
    pKa and pH."""

    ROUNDING_DIGITS = 5

    @staticmethod
    def calculate_acid_charge(pKa: float, pH: float) -> float:
        return round(-1.0 / (1.0 + 10 ** (pKa - pH)), ChargeCalculator.ROUNDING_DIGITS)

    @staticmethod
    def calculate_base_charge(pKa: float, pH: float) -> float:
        return round(1.0 / (1.0 + 10 ** (pH - pKa)), ChargeCalculator.ROUNDING_DIGITS)


if __name__ == "__main__":
    # Test
    example = {
        "pka_acid1": 3,
        "pka_acid2": 5.5,
        "pka_base1": None,
        "pka_base2": None,
    }
    print(pIPredictor.predict_input(example))
