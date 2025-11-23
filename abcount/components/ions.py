from enum import Enum

from abcount.components.classifier import ABClassData
from abcount.model.classifier import AcidType
from abcount.model.classifier import BaseType
from abcount.model.ions import IonDefinition


class IonMatcher:
    """Matches AB data to the ion and species definitions."""

    def match_class_data(self, class_data: ABClassData) -> IonDefinition:
        for ion_definition in IonRules:
            ion_obj = ion_definition.value
            if class_data == ion_obj.class_data:
                return ion_obj
        raise ValueError(f"Could not match any ion definition - (got {class_data})")


class IonRules(Enum):
    I1 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=AcidType.STRONG,
            acid_2_class=None,
            base_1_class=BaseType.STRONG,
            base_2_class=None,
        ),
        major_species_ph74_class="zwitterion",
        ion_class="zwitterion",
        explanation="zwitterion",
    )

    I2 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=AcidType.STRONG,
            acid_2_class=None,
            base_1_class=BaseType.WEAK,
            base_2_class=None,
        ),
        major_species_ph74_class="zwitterion",
        ion_class="zwitterion",
        explanation="zwitterion - weak base",
    )

    I3 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=AcidType.STRONG,
            acid_2_class=AcidType.STRONG,
            base_1_class=None,
            base_2_class=None,
        ),
        major_species_ph74_class="dianion",
        ion_class="diacid",
        explanation="diacid",
    )

    I4 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=AcidType.STRONG,
            acid_2_class=AcidType.WEAK,
            base_1_class=None,
            base_2_class=None,
        ),
        major_species_ph74_class="dianion",
        ion_class="diacid",
        explanation="diacid - strong + weak",
    )

    I5 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=AcidType.STRONG,
            acid_2_class=None,
            base_1_class=None,
            base_2_class=None,
        ),
        major_species_ph74_class="anion",
        ion_class="acid",
        explanation="acid",
    )

    I6 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=AcidType.WEAK,
            acid_2_class=None,
            base_1_class=None,
            base_2_class=None,
        ),
        major_species_ph74_class="anion",
        ion_class="weak acid",
        explanation="weak acid",
    )

    I7 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=AcidType.WEAK,
            acid_2_class=AcidType.WEAK,
            base_1_class=None,
            base_2_class=None,
        ),
        major_species_ph74_class="dianion",
        ion_class="weak diacid",
        explanation="weak diacid",
    )

    I8 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=None,
            acid_2_class=None,
            base_1_class=BaseType.STRONG,
            base_2_class=None,
        ),
        major_species_ph74_class="cation",
        ion_class="base",
        explanation="base",
    )

    I9 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=None,
            acid_2_class=None,
            base_1_class=BaseType.STRONG,
            base_2_class=BaseType.STRONG,
        ),
        major_species_ph74_class="dication",
        ion_class="dibase",
        explanation="dibase",
    )

    I10 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=None,
            acid_2_class=None,
            base_1_class=BaseType.WEAK,
            base_2_class=None,
        ),
        major_species_ph74_class="cation",
        ion_class="weak base",
        explanation="weak base",
    )

    I11 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=None,
            acid_2_class=None,
            base_1_class=BaseType.WEAK,
            base_2_class=BaseType.WEAK,
        ),
        major_species_ph74_class="dication",
        ion_class="weak dibase",
        explanation="weak dibase",
    )

    I12 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=None,
            acid_2_class=None,
            base_1_class=BaseType.STRONG,
            base_2_class=BaseType.WEAK,
        ),
        major_species_ph74_class="dication",
        ion_class="dibase",
        explanation="dibase - strong + weak",
    )

    I13 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=AcidType.WEAK,
            acid_2_class=None,
            base_1_class=BaseType.WEAK,
            base_2_class=None,
        ),
        major_species_ph74_class="zwitterion",
        ion_class="weak zwitterion",
        explanation="weak zwitterion",
    )

    I14 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=AcidType.WEAK,
            acid_2_class=None,
            base_1_class=BaseType.STRONG,
            base_2_class=None,
        ),
        major_species_ph74_class="zwitterion",
        ion_class="zwitterion",
        explanation="zwitterion - weak acid",
    )

    I15 = IonDefinition(
        class_data=ABClassData(
            acid_1_class=None, acid_2_class=None, base_1_class=None, base_2_class=None
        ),
        major_species_ph74_class="neutral",
        ion_class="neutral",
        explanation="neutral",
    )
