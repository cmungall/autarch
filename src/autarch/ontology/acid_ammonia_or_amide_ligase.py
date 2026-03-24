"""acid-ammonia (or amide) ligase.

Catalysis of the ligation of an acid to ammonia (NH3) or an amide via a
carbon-nitrogen bond, with the concomitant hydrolysis of the diphosphate bond
in ATP or a similar triphosphate.

Operational note: GO descendants placed under GO:0016880 include a broader set
of triphosphate-coupled acid-to-nitrogen ligations than the literal label
suggests, including small-molecule amines and some polymer-linked amide
formations. This classifier therefore models the shared chemistry: carboxylate
activation coupled to C-N bond formation with a non-carboxylate nitrogen
acceptor, while excluding amino-acid acceptors handled by EC 6.3.2. This is a
deliberately precision-biased approximation, so acid-bearing polyamine
acceptors that blur into peptide-ligase chemistry are left out rather than
absorbed with ad hoc exceptions.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from rdkit import Chem

from autarch.datamodel import ClassificationResult, Participant, PolymerType, Reaction
from autarch.moiety import is_carboxylic_acid
from autarch.molecules import (
    CHEBI_ADP,
    CHEBI_AMP,
    CHEBI_ATP,
    CHEBI_CDP,
    CHEBI_CMP,
    CHEBI_COA,
    CHEBI_CTP,
    CHEBI_DIPHOSPHATE,
    CHEBI_GDP,
    CHEBI_GTP,
    CHEBI_H2O,
    CHEBI_H_PLUS,
    CHEBI_NH3,
    CHEBI_NH4,
    CHEBI_O2,
    CHEBI_PHOSPHATE,
    CHEBI_UDP,
    CHEBI_UMP,
    CHEBI_UTP,
)
from autarch.ontology.ligase import Ligase

CHEBI_POLYPHOSPHATE = "CHEBI:16838"


class AcidAmmoniaOrAmideLigase(Ligase):
    """acid-ammonia (or amide) ligase.

    Catalysis of the ligation of an acid to ammonia (NH3) or an amide via a
    carbon-nitrogen bond, with the concomitant hydrolysis of the diphosphate
    bond in ATP or a similar triphosphate.
    """

    GO_ID = "GO:0016880"
    EC_NUMBER_PREFIX = "6.3.1.-"

    NTP_REACTANTS = {CHEBI_ATP, CHEBI_GTP, CHEBI_UTP, CHEBI_CTP}
    SPENT_NUCLEOTIDES = {
        CHEBI_ADP,
        CHEBI_AMP,
        CHEBI_GDP,
        CHEBI_UDP,
        CHEBI_UMP,
        CHEBI_CDP,
        CHEBI_CMP,
    }
    LEFT_SPECTATOR_CHEBIS = NTP_REACTANTS | {CHEBI_H2O, CHEBI_H_PLUS}
    RIGHT_SPECTATOR_CHEBIS = SPENT_NUCLEOTIDES | {
        CHEBI_PHOSPHATE,
        CHEBI_DIPHOSPHATE,
        CHEBI_POLYPHOSPHATE,
        CHEBI_H2O,
        CHEBI_H_PLUS,
    }
    AMIDE_PATTERN = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    FREE_AMINE_PATTERN = Chem.MolFromSmarts("[N;X3,X4+;H1,H2,H3+;!$(N[C,S,P]=O)]")
    AMIDE_NH_PATTERN = Chem.MolFromSmarts("[NX3;H1,H2][CX3](=[OX1])")
    POLYMER_NUCLEOPHILE_TYPES = {
        PolymerType.PROTEIN,
        PolymerType.PEPTIDE,
        PolymerType.PEPTIDOGLYCAN,
    }

    def check_membership_impl(self, reaction: Reaction) -> ClassificationResult:
        """Check for triphosphate-coupled acid-to-nitrogen ligation."""
        if reaction.is_transport_reaction():
            return ClassificationResult(
                is_member=False,
                explanation="Transport reaction - not acid-ammonia (or amide) ligation",
            )

        forward_result = self._check_direction(
            reaction.left_participants,
            reaction.right_participants,
        )
        if forward_result.is_member:
            return forward_result

        reverse_result = self._check_direction(
            reaction.right_participants,
            reaction.left_participants,
        )
        if reverse_result.is_member:
            return ClassificationResult(
                is_member=True,
                explanation=f"{reverse_result.explanation} (reverse reaction orientation)",
            )

        return forward_result

    def _check_direction(
        self,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> ClassificationResult:
        """Evaluate one reaction orientation."""
        if not self._has_triphosphate_coupling(left_participants, right_participants):
            return ClassificationResult(
                is_member=False,
                explanation="No ATP-or-similar triphosphate coupling detected",
            )

        if any(
            participant.chebi_id in {CHEBI_COA, CHEBI_O2}
            for participant in left_participants + right_participants
        ):
            return ClassificationResult(
                is_member=False,
                explanation="CoA or O2 chemistry indicates a different ligase branch",
            )

        left_reactive = [
            participant
            for participant in left_participants
            if participant.chebi_id not in self.LEFT_SPECTATOR_CHEBIS
        ]
        right_reactive = [
            participant
            for participant in right_participants
            if participant.chebi_id not in self.RIGHT_SPECTATOR_CHEBIS
        ]

        if len(left_reactive) < 2:
            return ClassificationResult(
                is_member=False,
                explanation="Requires both an acid donor and a nitrogen acceptor",
            )

        acid_donors = [
            participant for participant in left_reactive if self._is_carboxylic_acid_like(participant)
        ]
        if not acid_donors:
            return ClassificationResult(
                is_member=False,
                explanation="No carboxylic-acid donor detected",
            )

        if acid_donors and all(self._is_bicarbonate_like(participant) for participant in acid_donors):
            return ClassificationResult(
                is_member=False,
                explanation="Inorganic carbon fixation/carbamoylation chemistry is outside the acid-ammonia ligase branch",
            )

        nitrogen_acceptors = [
            participant
            for participant in left_reactive
            if participant not in acid_donors and self._is_nitrogen_acceptor(participant)
        ]
        if not nitrogen_acceptors:
            return ClassificationResult(
                is_member=False,
                explanation="No non-carboxylate nitrogen acceptor detected",
            )

        left_amide_bonds = self._amide_bond_count(left_reactive)
        right_amide_bonds = self._amide_bond_count(right_reactive)
        if right_amide_bonds > left_amide_bonds:
            return ClassificationResult(
                is_member=True,
                explanation=(
                    "Acid-ammonia (or amide) ligase: triphosphate-coupled amide formation "
                    "between a carboxylic-acid donor and a non-carboxylate nitrogen acceptor"
                ),
            )

        if self._has_polymer_ligation_signature(
            acid_donors=acid_donors,
            nitrogen_acceptors=nitrogen_acceptors,
            left_reactive=left_reactive,
            right_reactive=right_reactive,
        ):
            return ClassificationResult(
                is_member=True,
                explanation=(
                    "Acid-ammonia (or amide) ligase: triphosphate-coupled polymer-linked "
                    "acid-to-nitrogen ligation"
                ),
            )

        return ClassificationResult(
            is_member=False,
            explanation="No new amide linkage is formed in the principal products",
        )

    @classmethod
    def _has_triphosphate_coupling(
        cls,
        left_participants: list[Participant],
        right_participants: list[Participant],
    ) -> bool:
        """Check for ATP-or-similar triphosphate consumption."""
        has_ntp = any(participant.chebi_id in cls.NTP_REACTANTS for participant in left_participants)
        has_spent_nucleotide = any(
            participant.chebi_id in cls.SPENT_NUCLEOTIDES for participant in right_participants
        )
        has_phosphate_product = any(
            participant.chebi_id in {CHEBI_PHOSPHATE, CHEBI_DIPHOSPHATE, CHEBI_POLYPHOSPHATE}
            for participant in right_participants
        )
        return has_ntp and has_spent_nucleotide and has_phosphate_product

    @classmethod
    def _has_polymer_ligation_signature(
        cls,
        acid_donors: list[Participant],
        nitrogen_acceptors: list[Participant],
        left_reactive: list[Participant],
        right_reactive: list[Participant],
    ) -> bool:
        """Fallback for polymer ligations whose generic participants lack full structure."""
        if not acid_donors or not nitrogen_acceptors:
            return False
        left_polymers = [p for p in left_reactive if p.polymer_type in cls.POLYMER_NUCLEOPHILE_TYPES]
        right_polymers = [p for p in right_reactive if p.polymer_type in cls.POLYMER_NUCLEOPHILE_TYPES]
        if not left_polymers or not right_polymers:
            return False
        if not any(p.polymer_type in cls.POLYMER_NUCLEOPHILE_TYPES for p in nitrogen_acceptors):
            return False
        return any(
            left_polymer.polymer_type == right_polymer.polymer_type
            and not left_polymer.is_same_molecule(right_polymer)
            for left_polymer in left_polymers
            for right_polymer in right_polymers
        )

    @classmethod
    def _amide_bond_count(cls, participants: list[Participant]) -> int:
        """Count amide bonds across participants using cached ChEBI SMILES when needed."""
        total = 0
        for participant in participants:
            mol = cls._resolved_mol(participant)
            if mol is None or cls.AMIDE_PATTERN is None:
                continue
            total += len(mol.GetSubstructMatches(cls.AMIDE_PATTERN)) * max(participant.count, 1)
        return total

    @classmethod
    def _is_carboxylic_acid_like(cls, participant: Participant) -> bool:
        """Check whether the participant contains a carboxylic acid group."""
        smiles = cls._resolved_smiles(participant)
        return bool(smiles and is_carboxylic_acid(smiles, participant.chebi_id))

    @classmethod
    def _is_nitrogen_acceptor(cls, participant: Participant) -> bool:
        """Check for ammonia, amines, amides, or polymeric nitrogen acceptors."""
        if cls._is_free_ammonia(participant):
            return True
        if participant.polymer_type in cls.POLYMER_NUCLEOPHILE_TYPES:
            return True
        mol = cls._resolved_mol(participant)
        if mol is None:
            return False
        is_carboxylate = cls._is_carboxylic_acid_like(participant)
        if cls.FREE_AMINE_PATTERN and mol.HasSubstructMatch(cls.FREE_AMINE_PATTERN):
            return not (is_carboxylate and cls._is_amino_acid_like(participant))
        if cls.AMIDE_NH_PATTERN and mol.HasSubstructMatch(cls.AMIDE_NH_PATTERN):
            if is_carboxylate and cls._is_amino_acid_like(participant):
                return False
            return cls._amide_bond_count([participant]) <= 2
        return False

    @classmethod
    def _is_free_ammonia(cls, participant: Participant) -> bool:
        """Check for free ammonia or ammonium."""
        if participant.chebi_id in {CHEBI_NH3, CHEBI_NH4}:
            return True
        mol = cls._resolved_mol(participant)
        if mol is None:
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        return carbon_count == 0 and nitrogen_count == 1

    @classmethod
    def _is_amino_acid_like(cls, participant: Participant) -> bool:
        """Check for a free amino-acid-like substrate that belongs in EC 6.3.2."""
        mol = cls._resolved_mol(participant)
        if mol is None or cls.FREE_AMINE_PATTERN is None:
            return False
        if not cls._is_carboxylic_acid_like(participant):
            return False
        if not mol.HasSubstructMatch(cls.FREE_AMINE_PATTERN):
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        return carbon_count <= 12 and cls._amide_bond_count([participant]) <= 1

    @classmethod
    def _is_bicarbonate_like(cls, participant: Participant) -> bool:
        """Check for inorganic carbon donors such as bicarbonate/carbonate."""
        mol = cls._resolved_mol(participant)
        if mol is None:
            return False
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
        return carbon_count == 1 and oxygen_count >= 3 and nitrogen_count == 0 and phosphorus_count == 0

    @classmethod
    def _resolved_mol(cls, participant: Participant) -> Chem.Mol | None:
        """Get an RDKit molecule, resolving missing SMILES from the ChEBI cache."""
        smiles = cls._resolved_smiles(participant)
        if not smiles:
            return None
        return Chem.MolFromSmiles(smiles)

    @classmethod
    def _resolved_smiles(cls, participant: Participant) -> str | None:
        """Resolve participant SMILES directly or via the local ChEBI cache."""
        if participant.smiles:
            return participant.smiles
        if participant.chebi_id:
            return cls._load_chebi_smiles().get(participant.chebi_id)
        return None

    @staticmethod
    @lru_cache(maxsize=1)
    def _load_chebi_smiles() -> dict[str, str]:
        """Load the cached ChEBI SMILES map once."""
        cache_path = Path(__file__).resolve().parents[3] / "cache" / "chebi_smiles.json"
        if not cache_path.exists():
            return {}
        with cache_path.open() as handle:
            data = json.load(handle)
        return data if isinstance(data, dict) else {}
