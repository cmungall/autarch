"""catalytic activity, acting on a protein."""

from autarch.ontology.go_aggregate import ExplicitGoAggregate
from autarch.ontology.aminopeptidase import Aminopeptidase
from autarch.ontology.biotin_biotin_carboxyl_carrier_protein_ligase import BiotinBiotinCarboxylCarrierProteinLigase
from autarch.ontology.peptidase import Peptidase
from autarch.ontology.peptide_o_fucosyltransferase import PeptideOFucosyltransferase
from autarch.ontology.peptidyl_prolyl_cis_trans_isomerase import PeptidylProlylCisTransIsomerase
from autarch.ontology.protein_deglycase import ProteinDeglycase
from autarch.ontology.sumo_transferase import SUMOTransferase
from autarch.ontology.ubiquitin_protein_ligase import UbiquitinProteinLigase


class CatalyticActivityActingOnAProtein(ExplicitGoAggregate):
    """catalytic activity, acting on a protein."""

    GO_ID = "GO:0140096"
    CONCEPT_PHRASE = "catalytic activity, acting on a protein"
    CHILD_CLASSES = (Aminopeptidase, Peptidase, ProteinDeglycase, SUMOTransferase, UbiquitinProteinLigase, BiotinBiotinCarboxylCarrierProteinLigase, PeptidylProlylCisTransIsomerase, PeptideOFucosyltransferase,)
