"""Amino acid dehydrogenase reaction classification.

EC 1.4.1: Oxidoreductases acting on CH-NH2 group with NAD(P)+ as acceptor.
These enzymes catalyze oxidative deamination of amino acids to 2-oxoacids.

Pattern: amino acid + NAD(P)+ + H2O = 2-oxoacid + NH4+ + NAD(P)H + H+
"""

from typing import ClassVar, Optional

from autarch.datamodel import ClassificationResult, Reaction
from autarch.ontology.oxidoreductase import Oxidoreductase
from autarch.molecules import (
    CHEBI_NAD_PLUS,
    CHEBI_NADH,
    CHEBI_NADP_PLUS,
    CHEBI_NADPH,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
    CHEBI_H2O2,
)


class OxidoreductaseActingOnTheCHNH2GroupOfDonorsNADOrNADPAsAcceptor(Oxidoreductase):
    """oxidoreductase acting on the CH-NH2 group of donors NAD or NADP as acceptor

    EC 1.4.1.x includes:
    - Glutamate dehydrogenase: L-glutamate + NAD(P)+ + H2O → 2-oxoglutarate + NH4+ + NAD(P)H
    - Alanine dehydrogenase: L-alanine + NAD+ + H2O → pyruvate + NH4+ + NADH
    - Leucine dehydrogenase: L-leucine + NAD+ + H2O → 4-methyl-2-oxopentanoate + NH4+ + NADH
    - Valine dehydrogenase: L-valine + NAD(P)+ + H2O → 3-methyl-2-oxobutanoate + NH4+ + NAD(P)H
    - Aspartate dehydrogenase: L-aspartate + NAD(P)+ + H2O → oxaloacetate + NH4+ + NAD(P)H

    Key distinction from EC 1.4.3 (amino acid oxidases):
    - EC 1.4.1 uses NAD(P)+ as electron acceptor
    - EC 1.4.3 uses O2 as electron acceptor (produces H2O2)

    Examples:
    - L-glutamate + NADP+ + H2O = 2-oxoglutarate + NH4+ + NADPH + H+
    - L-alanine + NAD+ + H2O = pyruvate + NH4+ + NADH + H+
    """

    GO_ID = "GO:0016639"  # oxidoreductase activity, acting on CH-NH2 with NAD(P)+
    EC_NUMBER_PREFIX: ClassVar[Optional[str]] = "1.4.1.-"

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check if reaction is an amino acid dehydrogenase.

        Strategy:
        1. Must use NAD(P)+ as substrate and produce NAD(P)H
        2. Must release NH4+/NH3 (deamination)
        3. Must NOT use O2 (that would be amino acid oxidase, EC 1.4.3)
        4. Look for characteristic amino acid substrates or 2-oxoacid products
        """
        # Use RHEA label for pattern matching (not p.name which has ChEBI generic names)
        label_lower = reaction.label.lower() if reaction.label else ""
        label_parts = label_lower.split("=")
        substrate_str = label_parts[0] if label_parts else ""
        product_str = label_parts[1] if len(label_parts) > 1 else ""

        # Must have NAD(P)+ as substrate
        has_nad_plus = any(
            p.chebi_id in {CHEBI_NAD_PLUS, CHEBI_NADP_PLUS}
            for p in reaction.left_participants
        )

        if not has_nad_plus:
            return ClassificationResult(
                is_member=False,
                explanation="No NAD(P)+ substrate"
            )

        # Must produce NAD(P)H
        has_nadh = any(
            p.chebi_id in {CHEBI_NADH, CHEBI_NADPH}
            for p in reaction.right_participants
        )

        if not has_nadh:
            return ClassificationResult(
                is_member=False,
                explanation="No NAD(P)H product"
            )

        # Must release ammonia/ammonium (deamination)
        has_ammonia_product = any(
            p.chebi_id in {CHEBI_NH3, CHEBI_NH4}
            for p in reaction.right_participants
        ) or any(x in product_str for x in ["ammonia", "ammonium", "nh4", "nh3"])

        if not has_ammonia_product:
            return ClassificationResult(
                is_member=False,
                explanation="No ammonia/ammonium release - not deamination"
            )

        # Must NOT use O2 (that would be amino acid oxidase, EC 1.4.3)
        has_o2 = any(
            p.chebi_id == CHEBI_O2
            for p in reaction.left_participants
        )

        if has_o2:
            return ClassificationResult(
                is_member=False,
                explanation="Uses O2 - amino acid oxidase (EC 1.4.3), not dehydrogenase"
            )

        # Must NOT produce H2O2 (another indicator of EC 1.4.3)
        has_h2o2 = any(
            p.chebi_id == CHEBI_H2O2
            for p in reaction.right_participants
        )

        if has_h2o2:
            return ClassificationResult(
                is_member=False,
                explanation="Produces H2O2 - amino acid oxidase (EC 1.4.3)"
            )

        # Characteristic amino acid substrates
        amino_acid_patterns = [
            "glutamate", "glutamic", "aspartate", "aspartic",
            "alanine", "valine", "leucine", "isoleucine",
            "lysine", "ornithine", "arginine",
            "tryptophan", "phenylalanine", "tyrosine",
            "serine", "threonine", "glycine",
            "proline", "histidine", "methionine",
            "cysteine", "diaminopimelate", "alpha-amino",
        ]

        has_amino_acid = any(aa in substrate_str for aa in amino_acid_patterns)

        # Characteristic 2-oxoacid products
        oxoacid_patterns = [
            "oxoglutarate", "pyruvate", "oxaloacetate",
            "oxobutanoate", "oxopentanoate", "oxoheptane",
            "oxovalerate", "oxocaproate", "oxohexan",
            "2-oxo", "α-oxo", "alpha-oxo", "keto",
            "oxoisocaproate", "oxoisovalerate", "indole-3-pyruvate",
            "piperideine-6-carboxylate",  # from lysine dehydrogenase
        ]

        has_oxoacid = any(oa in product_str for oa in oxoacid_patterns)

        # Must have either amino acid substrate OR oxoacid product
        if has_amino_acid or has_oxoacid:
            # Derive substrate name from label
            substrate_name = next(
                (aa for aa in amino_acid_patterns if aa in substrate_str),
                "amino acid"
            )

            return ClassificationResult(
                is_member=True,
                explanation=f"Amino acid dehydrogenase: {substrate_name} + NAD(P)+ → 2-oxoacid + NH4+"
            )

        return ClassificationResult(
            is_member=False,
            explanation="No amino acid dehydrogenase pattern"
        )
