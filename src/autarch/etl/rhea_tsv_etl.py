"""RHEA TSV ETL - Clean, simple alternative to SPARQL approach.

This module provides a straightforward TSV-based ETL for RHEA reactions,
eliminating the complexity and duplication bugs of the SPARQL approach.
Uses authoritative TSV files directly from RHEA's FTP site.

Supports polymer reactions with (n) stoichiometry notation by parsing
human-readable labels from rhea-reactions.txt.gz.
"""

import gzip
import logging
import pandas as pd
import re
import requests
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from urllib.parse import urljoin

from autarch.datamodel import Reaction, Participant, RheaTerm, PolymerType
from autarch.etl.label_parser import (
    parse_reaction_label,
    has_polymer_notation,
    ParsedParticipant,
)
from autarch.etl.chebi_smiles import canonicalize_smiles
from autarch.etl.chebi_rdkit_cache import get_inchi
from autarch.etl.chebi_normalization import (
    EquationChEBILookup,
    normalize_chebi_by_name,
)
from rdkit import Chem

logger = logging.getLogger(__name__)

# Lazy-loaded ChEBI lookup cache
_SMILES_TO_CHEBI: Optional[Dict[str, str]] = None
_CHEBI_TO_NAME: Optional[Dict[str, str]] = None


def _load_chebi_lookup() -> Tuple[Dict[str, str], Dict[str, str]]:
    """Load or return cached SMILES→ChEBI lookup.

    Returns:
        Tuple of (smiles_to_chebi, chebi_to_name) dictionaries
    """
    global _SMILES_TO_CHEBI, _CHEBI_TO_NAME

    if _SMILES_TO_CHEBI is None:
        import json
        from pathlib import Path

        cache_dir = Path("cache")
        smiles_cache = cache_dir / "smiles_to_chebi.json"
        names_cache = cache_dir / "chebi_names.json"

        if smiles_cache.exists():
            with open(smiles_cache) as f:
                _SMILES_TO_CHEBI = json.load(f)
            logger.info(f"Loaded {len(_SMILES_TO_CHEBI)} SMILES→ChEBI mappings")
        else:
            logger.warning(
                "SMILES→ChEBI cache not found. Run 'autarch cache-smiles-lookup' to build it."
            )
            _SMILES_TO_CHEBI = {}

        if names_cache.exists():
            with open(names_cache) as f:
                _CHEBI_TO_NAME = json.load(f)
            logger.info(f"Loaded {len(_CHEBI_TO_NAME)} ChEBI names")
        else:
            _CHEBI_TO_NAME = {}

    assert _SMILES_TO_CHEBI is not None
    assert _CHEBI_TO_NAME is not None
    return _SMILES_TO_CHEBI, _CHEBI_TO_NAME


def lookup_chebi_for_smiles(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Look up ChEBI ID and name for a SMILES string.

    Args:
        smiles: SMILES string (will be canonicalized)

    Returns:
        Tuple of (chebi_id, name) or (None, None) if not found
    """
    smiles_to_chebi, chebi_to_name = _load_chebi_lookup()

    canonical = canonicalize_smiles(smiles)
    if canonical is None:
        return None, None

    chebi_id = smiles_to_chebi.get(canonical)
    if chebi_id is None:
        return None, None

    name = chebi_to_name.get(chebi_id)
    return chebi_id, name

# RHEA FTP base URL for TSV files
RHEA_FTP_BASE = "https://ftp.expasy.org/databases/rhea/tsv/"
RHEA_TXT_BASE = "https://ftp.expasy.org/databases/rhea/txt/"

# Core TSV files we need
RHEA_TSV_FILES = {
    'reactions': 'rhea-reaction-smiles.tsv',        # Reaction SMILES (main file)
    'participants': 'chebiId_name.tsv',             # ChEBI ID to name mapping
    'directions': 'rhea-directions.tsv',            # Reaction direction info
    'ec_mapping': 'rhea2ec.tsv',                    # EC number mappings
    'chebi_smiles': 'rhea-chebi-smiles.tsv',        # ChEBI SMILES
}

# Additional text files for labels (needed for polymer stoichiometry)
RHEA_TXT_FILES = {
    'reactions_txt': 'rhea-reactions.txt.gz',       # Human-readable labels with (n) notation
}


@dataclass
class RheaTSVETL:
    """Simple, clean RHEA ETL using TSV files."""
    
    cache_dir: Path = Path("cache/rhea_tsv")
    force_download: bool = False
    
    def __post_init__(self):
        """Ensure cache directory exists."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
    def download_tsv_file(self, filename: str) -> Path:
        """Download a TSV file from RHEA FTP site.
        
        Args:
            filename: Name of TSV file to download
            
        Returns:
            Path to downloaded file
        """
        local_path = self.cache_dir / filename
        
        if local_path.exists() and not self.force_download:
            logger.info(f"Using cached file: {local_path}")
            return local_path
            
        url = urljoin(RHEA_FTP_BASE, filename)
        logger.info(f"Downloading {filename} from {url}")
        
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            with open(local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
                    
            logger.info(f"Downloaded {filename} to {local_path}")
            return local_path
            
        except Exception as e:
            logger.error(f"Failed to download {filename}: {e}")
            raise
            
    def download_all_files(self) -> Dict[str, Path]:
        """Download all required TSV files.

        Returns:
            Dictionary mapping file types to local paths
        """
        paths = {}
        for file_type, filename in RHEA_TSV_FILES.items():
            try:
                paths[file_type] = self.download_tsv_file(filename)
            except Exception as e:
                logger.warning(f"Could not download {file_type} ({filename}): {e}")

        return paths

    def download_reactions_txt(self) -> Path:
        """Download rhea-reactions.txt.gz for human-readable labels.

        This file contains DEFINITION lines with (n) polymer notation
        that we need for proper stoichiometry parsing.

        Returns:
            Path to downloaded file
        """
        filename = RHEA_TXT_FILES['reactions_txt']
        local_path = self.cache_dir / filename

        if local_path.exists() and not self.force_download:
            logger.info(f"Using cached file: {local_path}")
            return local_path

        url = urljoin(RHEA_TXT_BASE, filename)
        logger.info(f"Downloading {filename} from {url}")

        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(local_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        logger.info(f"Downloaded {filename} to {local_path}")
        return local_path

    def load_reaction_labels(self) -> Dict[str, str]:
        """Load human-readable reaction labels from rhea-reactions.txt.gz.

        Parses the KEGG-style format to extract RHEA ID -> DEFINITION mapping.
        The DEFINITION contains the human-readable equation with (n) notation.

        Returns:
            Dictionary mapping RHEA IDs (without prefix) to labels

        Examples:
            >>> etl = RheaTSVETL()
            >>> labels = etl.load_reaction_labels()
            >>> labels.get('10256')
            '[(1->4)-6-phospho-alpha-D-glucosyl](n) + n ATP + n H2O = ...'
        """
        txt_path = self.download_reactions_txt()

        labels = {}
        current_rhea_id = None

        with gzip.open(txt_path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.rstrip()

                if line.startswith('ENTRY'):
                    # Parse RHEA ID: "ENTRY       RHEA:10000"
                    match = re.match(r'ENTRY\s+RHEA:(\d+)', line)
                    if match:
                        current_rhea_id = match.group(1)

                elif line.startswith('DEFINITION') and current_rhea_id:
                    # Parse definition: "DEFINITION  pentanamide + H2O = ..."
                    definition = line[12:].strip()  # Skip "DEFINITION  "
                    labels[current_rhea_id] = definition

                elif line.startswith('///'):
                    current_rhea_id = None

        logger.info(f"Loaded {len(labels)} reaction labels")
        return labels

    def get_polymer_labels(self) -> Dict[str, str]:
        """Get labels for reactions that have polymer participants.

        This includes:
        1. Reactions with (n) notation for variable stoichiometry
        2. Reactions with polymer-type molecules (tRNA, protein, DNA, RNA)
           even without (n) notation

        Returns:
            Dictionary mapping RHEA IDs to labels, for polymer reactions
        """
        all_labels = self.load_reaction_labels()

        # Keywords indicating polymer-type participants
        polymer_type_keywords = [
            'trna', '[protein]', 'protein]', '-protein]',
            '[dna', 'dna]', '[rna', 'rna(',
            'mrna', 'rrna', 'snrna', 'lncrna',
            'polysaccharide', 'glycogen', 'starch',
            'chitin', 'cellulose', 'peptidoglycan',
        ]

        polymer_labels = {}
        for rhea_id, label in all_labels.items():
            label_lower = label.lower()
            # Include if has (n) notation
            if has_polymer_notation(label):
                polymer_labels[rhea_id] = label
            # Include if has polymer-type keywords
            elif any(kw in label_lower for kw in polymer_type_keywords):
                polymer_labels[rhea_id] = label

        logger.info(f"Found {len(polymer_labels)} reactions with polymer participants")
        return polymer_labels

    def load_equation_chebi_lookup(self) -> EquationChEBILookup:
        """Load ChEBI IDs from EQUATION field in rhea-reactions.txt.gz.

        The EQUATION field contains authoritative ChEBI IDs that can be used
        to fill in missing IDs when SMILES lookup fails.

        Returns:
            EquationChEBILookup instance with loaded data

        Examples:
            >>> etl = RheaTSVETL()  # doctest: +SKIP
            >>> lookup = etl.load_equation_chebi_lookup()  # doctest: +SKIP
            >>> left, right = lookup.get_chebi_ids("10000")  # doctest: +SKIP
        """
        txt_path = self.download_reactions_txt()
        lookup = EquationChEBILookup()
        count = lookup.load_from_file(str(txt_path))
        logger.info(f"Loaded {count} EQUATION ChEBI mappings")
        return lookup

    def load_directions(self) -> Dict[str, str]:
        """Load rhea-directions.tsv to get proper master ID mappings.

        The directions file maps each reaction ID to its master/undefined ID.
        This is more reliable than the arithmetic formula (id // 4) * 4.

        Returns:
            Dictionary mapping any RHEA ID to its master ID (both as strings)

        Examples:
            >>> etl = RheaTSVETL()
            >>> dirs = etl.load_directions()  # doctest: +SKIP
            >>> dirs.get("12818")  # doctest: +SKIP
            '12817'
        """
        directions_path = self.download_tsv_file(RHEA_TSV_FILES['directions'])

        # Maps any reaction ID to its master ID
        id_to_master: Dict[str, str] = {}

        with open(directions_path, 'r') as f:
            # Skip header: RHEA_ID_MASTER, RHEA_ID_LR, RHEA_ID_RL, RHEA_ID_BI
            f.readline()

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    master_id = parts[0]  # RHEA_ID_MASTER (undefined direction)
                    lr_id = parts[1]      # RHEA_ID_LR
                    rl_id = parts[2]      # RHEA_ID_RL
                    bi_id = parts[3]      # RHEA_ID_BI (bidirectional)

                    # Map all variants to the master ID
                    for variant_id in [master_id, lr_id, rl_id, bi_id]:
                        if variant_id:
                            id_to_master[variant_id] = master_id

        logger.info(f"Loaded direction mappings for {len(id_to_master)} reaction IDs")
        return id_to_master

    def parsed_to_participant(self, parsed: ParsedParticipant) -> Participant:
        """Convert a ParsedParticipant from label parsing to a Participant object.

        Args:
            parsed: ParsedParticipant from label_parser

        Returns:
            Participant object with polymer information populated
        """
        # Determine polymer type enum
        polymer_type_enum = None
        if parsed.polymer_type:
            try:
                polymer_type_enum = PolymerType(parsed.polymer_type)
            except ValueError:
                polymer_type_enum = PolymerType.OTHER

        return Participant(
            name=parsed.base_name or parsed.name,
            stoichiometry=parsed.stoichiometry,
            polymer_index=parsed.polymer_index,
            polymer_type=polymer_type_enum,
            linkage=parsed.linkage,
            monomer=parsed.monomer,
            location=parsed.location,
        )

    def reaction_from_label(self, rhea_id: str, label: str) -> Optional[Reaction]:
        """Create a Reaction object from a label string (no SMILES required).

        This enables representing polymer reactions that cannot be expressed
        as SMILES strings.

        Args:
            rhea_id: RHEA reaction ID
            label: Human-readable reaction equation

        Returns:
            Reaction object or None if parsing failed

        Examples:
            >>> etl = RheaTSVETL()
            >>> rxn = etl.reaction_from_label("10256", "RNA(n) + NTP = RNA(n+1) + PPi")
            >>> rxn is not None
            True
            >>> len(rxn.left_participants)
            2
        """
        left_parsed, right_parsed = parse_reaction_label(label)

        if not left_parsed and not right_parsed:
            logger.warning(f"Could not parse label for {rhea_id}: {label[:50]}...")
            return None

        left_participants = [self.parsed_to_participant(p) for p in left_parsed]
        right_participants = [self.parsed_to_participant(p) for p in right_parsed]

        return Reaction(
            left_participants=left_participants,
            right_participants=right_participants,
        )

    def load_polymer_reactions(self) -> Dict[str, Reaction]:
        """Load polymer reactions from labels (no SMILES required).

        These are reactions with (n) notation that cannot be represented
        as SMILES strings but can be parsed from human-readable labels.

        Returns:
            Dictionary mapping RHEA IDs to Reaction objects
        """
        polymer_labels = self.get_polymer_labels()

        # Get IDs that already have SMILES (to avoid duplicates)
        reactions_df = self.load_reaction_smiles()
        smiles_ids = set(reactions_df['rhea_id'].astype(str))

        # Only include polymer reactions that don't have SMILES
        label_only_ids = set(polymer_labels.keys()) - smiles_ids

        logger.info(f"Loading {len(label_only_ids)} polymer reactions from labels (no SMILES)")

        reactions = {}
        for rhea_id in label_only_ids:
            label = polymer_labels[rhea_id]
            reaction = self.reaction_from_label(rhea_id, label)
            if reaction:
                reactions[rhea_id] = reaction

        logger.info(f"Successfully loaded {len(reactions)} polymer reactions from labels")
        return reactions
        
    def load_reaction_smiles(self) -> pd.DataFrame:
        """Load reaction SMILES data.
        
        Returns:
            DataFrame with columns: rhea_id, reaction_smiles
        """
        reactions_path = self.download_tsv_file(RHEA_TSV_FILES['reactions'])
        
        logger.info(f"Loading reaction SMILES from {reactions_path}")
        
        # Load TSV - simple, clean data!
        df = pd.read_csv(
            reactions_path, 
            sep='\t', 
            header=None,
            names=['rhea_id', 'reaction_smiles'],
            dtype={'rhea_id': str, 'reaction_smiles': str}
        )
        
        logger.info(f"Loaded {len(df)} reaction SMILES")
        return df
        
    def load_chebi_names(self) -> Dict[str, str]:
        """Load ChEBI ID to name mappings.
        
        Returns:
            Dictionary mapping ChEBI IDs to names
        """
        try:
            names_path = self.download_tsv_file(RHEA_TSV_FILES['participants'])
            
            df = pd.read_csv(
                names_path,
                sep='\t', 
                header=None,
                names=['chebi_id', 'name'],
                dtype=str
            )
            
            # Convert to dictionary for fast lookup
            name_map = dict(zip(df['chebi_id'], df['name']))
            logger.info(f"Loaded {len(name_map)} ChEBI names")
            return name_map
            
        except Exception as e:
            logger.warning(f"Could not load ChEBI names: {e}")
            return {}
            
    def load_go_mappings(self, go_cache_path: Optional[Path] = None) -> Dict[int, List[str]]:
        """Load GO term mappings from GO cache.

        Returns mapping from RHEA ID to list of GO terms.
        GO is the source of truth for RHEA→GO mappings.

        Note: GO stores master RHEA IDs directly (not computed).
        RHEA master IDs are NOT always divisible by 4 - they're defined
        in the rhea-directions.tsv file.

        Args:
            go_cache_path: Path to go_terms.jsonl file

        Returns:
            Dictionary mapping RHEA IDs (as int) to list of GO term IDs
        """
        import json

        if go_cache_path is None:
            go_cache_path = Path("cache") / "go_terms.jsonl"

        if not go_cache_path.exists():
            logger.warning(f"GO cache not found at {go_cache_path}")
            return {}

        # Build RHEA ID -> GO terms mapping
        # GO stores master RHEA IDs directly - don't compute them
        rhea_to_go: Dict[int, List[str]] = {}

        with open(go_cache_path) as f:
            for line in f:
                term_data = json.loads(line)
                go_id = term_data.get("go_id")
                for rhea_id_str in term_data.get("rhea_ids", []):
                    if rhea_id_str.startswith("RHEA:"):
                        # Use the ID as-is from GO (GO stores master IDs)
                        numeric_id = int(rhea_id_str[5:])
                        if numeric_id not in rhea_to_go:
                            rhea_to_go[numeric_id] = []
                        if go_id not in rhea_to_go[numeric_id]:
                            rhea_to_go[numeric_id].append(go_id)

        logger.info(f"Loaded GO mappings for {len(rhea_to_go)} RHEA IDs from GO (source of truth)")
        return rhea_to_go

    def load_ec_mappings(self) -> Dict[str, List[str]]:
        """Load EC number mappings.

        RHEA uses a hierarchical ID structure:
        - Master ID (e.g., 10000) is always divisible by 4
        - Directional IDs: master+1 (LR), master+2 (RL), master+3 (BI)

        The EC file only contains master IDs, but SMILES reactions use
        directional IDs. This method maps both master and directional IDs
        to their EC numbers.

        Returns:
            Dictionary mapping RHEA IDs (both master and directional) to EC numbers
        """
        try:
            ec_path = self.download_tsv_file(RHEA_TSV_FILES['ec_mapping'])

            df = pd.read_csv(
                ec_path,
                sep='\t',
                dtype=str
            )

            # The file has columns: RHEA_ID, DIRECTION, MASTER_ID, ID
            # Group by MASTER_ID to get all EC numbers for a reaction
            master_ec_map = df.groupby('MASTER_ID')['ID'].apply(list).to_dict()

            # Create mappings for both master and directional IDs
            ec_map = {}
            for master_id, ec_numbers in master_ec_map.items():
                master_int = int(master_id)
                # Add mapping for master ID
                ec_map[master_id] = ec_numbers
                # Add mappings for directional IDs (master+1, +2, +3)
                for offset in [1, 2, 3]:
                    directional_id = str(master_int + offset)
                    ec_map[directional_id] = ec_numbers

            logger.info(f"Loaded EC mappings for {len(master_ec_map)} master reactions "
                       f"({len(ec_map)} total including directional variants)")
            return ec_map

        except Exception as e:
            logger.warning(f"Could not load EC mappings: {e}")
            return {}
            
    def smiles_to_participants(
        self,
        smiles: str,
        is_left_side: bool = True,
        chebi_names: Optional[Dict[str, str]] = None
    ) -> List[Participant]:
        """Convert SMILES side to Participant objects.

        Args:
            smiles: SMILES string for one side of reaction
            is_left_side: Whether this is left side (True) or right side (False)
            chebi_names: Optional ChEBI name lookup (legacy, now uses cache)

        Returns:
            List of Participant objects
        """
        if not smiles or smiles.strip() == '':
            return []

        participants = []

        # Split on '.' to get individual molecules
        molecule_smiles = smiles.split('.')

        for mol_smiles in molecule_smiles:
            mol_smiles = mol_smiles.strip()
            if not mol_smiles:
                continue

            # Handle polymer wildcards (*) - these represent protein/polymer parts
            if '*' in mol_smiles:
                name = self._get_polymer_name(mol_smiles)
                participant = Participant(
                    smiles=mol_smiles,
                    name=name,
                    count=1  # TSV already has correct stoichiometry
                )
            else:
                # Look up ChEBI ID and name from cache
                chebi_id, chebi_name = lookup_chebi_for_smiles(mol_smiles)

                # Use ChEBI name if available, otherwise guess from SMILES
                if chebi_name:
                    name = chebi_name
                else:
                    name = self._guess_molecule_name(mol_smiles, chebi_names)

                # Validate with RDKit
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol is None:
                    logger.debug(f"Could not parse SMILES: {mol_smiles}")

                # Get InChI for stereochemistry comparison (cached by ChEBI ID)
                inchi = None
                if chebi_id and mol_smiles:
                    inchi = get_inchi(chebi_id, mol_smiles)

                # Create participant with ChEBI ID, name, SMILES, and InChI
                participant = Participant(
                    chebi_id=chebi_id,
                    smiles=mol_smiles,
                    inchi=inchi,
                    name=name,
                    count=1
                )

            participants.append(participant)

        return participants
        
    def _get_polymer_name(self, smiles: str) -> str:
        """Generate name for polymer/protein from SMILES with wildcards.
        
        Args:
            smiles: SMILES string containing * wildcards
            
        Returns:
            Descriptive name for the polymer
        """
        # Common patterns in RHEA polymer SMILES
        if '*N[C@@H](' in smiles and 'C(*)=O' in smiles:
            if 'CS' in smiles:
                return "cysteine-containing protein"
            elif 'COP' in smiles:
                return "phosphorylated protein"  
            else:
                return "amino acid protein"
        elif '*OO' in smiles:
            return "peroxide-containing compound"
        elif '*[C@H]' in smiles:
            return "carbohydrate polymer"
        else:
            return "generic polymer"
            
    def _guess_molecule_name(
        self, 
        smiles: str, 
        chebi_names: Optional[Dict[str, str]] = None
    ) -> str:
        """Guess molecule name from SMILES.
        
        Args:
            smiles: SMILES string
            chebi_names: Optional ChEBI name lookup
            
        Returns:
            Best guess at molecule name
        """
        # Common small molecules in biochemistry
        common_molecules = {
            'O': 'water',
            '[H+]': 'proton', 
            '[H]O[H]': 'water',
            'O=O': 'oxygen',
            'O=C=O': 'carbon dioxide',
            '[NH4+]': 'ammonium',
            'O=P([O-])([O-])[O-]': 'phosphate',
        }
        
        if smiles in common_molecules:
            return common_molecules[smiles]
            
        # For complex molecules, use a truncated SMILES as name
        return f"compound_{smiles[:15]}"

    def normalize_participant_chebi(
        self,
        participant: Participant,
        equation_chebi_ids: Optional[List[str]] = None,
        position: int = 0,
    ) -> Participant:
        """Normalize ChEBI ID for a participant using multiple fallback strategies.

        Strategy order:
        1. Keep existing ChEBI ID if present
        2. Try EQUATION ChEBI ID by position
        3. Try name-based normalization

        Args:
            participant: Participant to normalize
            equation_chebi_ids: List of ChEBI IDs from EQUATION (in order)
            position: Position of this participant in the reaction side

        Returns:
            Participant with normalized ChEBI ID (may be same as input)

        Examples:
            >>> etl = RheaTSVETL()
            >>> p = Participant(name="H(+)", chebi_id=None)
            >>> normalized = etl.normalize_participant_chebi(p)
            >>> normalized.chebi_id
            'CHEBI:15378'
        """
        # Strategy 1: Already has ChEBI ID
        if participant.chebi_id:
            return participant

        # Strategy 2: Use EQUATION ChEBI ID by position
        if equation_chebi_ids and position < len(equation_chebi_ids):
            participant.chebi_id = equation_chebi_ids[position]
            return participant

        # Strategy 3: Name-based normalization
        if participant.name:
            chebi_id = normalize_chebi_by_name(participant.name)
            if chebi_id:
                participant.chebi_id = chebi_id
                return participant

        return participant

    def normalize_reaction_chebi(
        self,
        reaction: Reaction,
        rhea_id: str,
        equation_lookup: Optional[EquationChEBILookup] = None,
    ) -> Reaction:
        """Normalize ChEBI IDs for all participants in a reaction.

        Uses EQUATION data and name-based normalization to fill in
        missing ChEBI IDs.

        Args:
            reaction: Reaction to normalize
            rhea_id: RHEA ID (without prefix) for EQUATION lookup
            equation_lookup: Optional pre-loaded EQUATION lookup

        Returns:
            Reaction with normalized ChEBI IDs
        """
        # Get EQUATION ChEBI IDs if available
        left_chebi_ids: List[str] = []
        right_chebi_ids: List[str] = []
        if equation_lookup:
            left_chebi_ids, right_chebi_ids = equation_lookup.get_chebi_ids(rhea_id)

        # Normalize left participants
        for i, participant in enumerate(reaction.left_participants):
            self.normalize_participant_chebi(participant, left_chebi_ids, i)

        # Normalize right participants
        for i, participant in enumerate(reaction.right_participants):
            self.normalize_participant_chebi(participant, right_chebi_ids, i)

        return reaction

    def reaction_smiles_to_reaction(
        self, 
        rhea_id: str, 
        reaction_smiles: str,
        chebi_names: Optional[Dict[str, str]] = None,
        ec_numbers: Optional[List[str]] = None
    ) -> Optional[Reaction]:
        """Convert reaction SMILES to Reaction object.
        
        Args:
            rhea_id: RHEA reaction ID
            reaction_smiles: Reaction SMILES string  
            chebi_names: Optional ChEBI name lookup
            ec_numbers: Optional EC numbers for this reaction
            
        Returns:
            Reaction object or None if parsing failed
        """
        if '>>' not in reaction_smiles:
            logger.debug(f"Invalid reaction SMILES for {rhea_id}: {reaction_smiles}")
            return None
            
        try:
            # Split reaction SMILES into left and right sides
            left_smiles, right_smiles = reaction_smiles.split('>>')
            
            # Convert each side to participants
            left_participants = self.smiles_to_participants(
                left_smiles, is_left_side=True, chebi_names=chebi_names
            )
            right_participants = self.smiles_to_participants(
                right_smiles, is_left_side=False, chebi_names=chebi_names  
            )
            
            # Create reaction object - CLEAN DATA, NO DUPLICATION!
            reaction = Reaction(
                left_participants=left_participants,
                right_participants=right_participants
            )
            
            logger.debug(f"Converted {rhea_id}: {len(left_participants)} → {len(right_participants)} participants")
            return reaction
            
        except Exception as e:
            logger.warning(f"Failed to parse reaction {rhea_id}: {e}")
            return None
            
    def load_reactions_batch(
        self,
        rhea_ids: Optional[List[str]] = None,
        limit: Optional[int] = None,
        include_polymer_reactions: bool = True,
        go_cache_path: Optional[Path] = None
    ) -> Dict[str, Reaction]:
        """Load reactions from TSV files and optionally from labels.

        Args:
            rhea_ids: Optional list of specific RHEA IDs to load
            limit: Optional limit on number of reactions to load
            include_polymer_reactions: Include polymer reactions parsed from labels
                (reactions with (n) notation that have no SMILES)
            go_cache_path: Path to GO cache (for including GO-mapped reactions)

        Returns:
            Dictionary mapping RHEA IDs to Reaction objects
        """
        logger.info("Loading reactions from TSV files...")

        # Load core data
        reactions_df = self.load_reaction_smiles()
        chebi_names = self.load_chebi_names()
        ec_mappings = self.load_ec_mappings()

        # Load GO mappings to include GO-mapped reactions (GO is source of truth)
        go_mappings = self.load_go_mappings(go_cache_path)
        go_rhea_ids = {str(rhea_id) for rhea_id in go_mappings.keys()}
        # Also include directional variants for GO-mapped reactions
        for master_id in list(go_rhea_ids):
            for offset in [1, 2, 3]:
                go_rhea_ids.add(str(int(master_id) + offset))

        # Filter if specific IDs requested
        if rhea_ids:
            reactions_df = reactions_df[reactions_df['rhea_id'].isin(rhea_ids)]
            logger.info(f"Filtered to {len(reactions_df)} requested reactions")
        else:
            # Include reactions that have EC annotations OR GO mappings
            # GO is the source of truth for mappings, so we should include those
            ec_reactions = set(ec_mappings.keys())
            relevant_reactions = ec_reactions | go_rhea_ids

            reactions_with_annotations = reactions_df[
                reactions_df['rhea_id'].isin(relevant_reactions)
            ]

            if len(reactions_with_annotations) > 0:
                ec_only = len(ec_reactions - go_rhea_ids)
                go_only = len(go_rhea_ids - ec_reactions)
                both = len(ec_reactions & go_rhea_ids)
                logger.info(
                    f"Found {len(reactions_with_annotations)} reactions with annotations "
                    f"(EC only: {ec_only}, GO only: {go_only}, both: {both})"
                )
                reactions_df = reactions_with_annotations
            else:
                logger.warning("No reactions found with annotations - using all reactions")

        # Limit if requested
        if limit:
            reactions_df = reactions_df.head(limit)
            logger.info(f"Limited to {len(reactions_df)} reactions")

        # Convert to Reaction objects
        reactions = {}
        failed_count = 0

        for _, row in reactions_df.iterrows():
            rhea_id = str(row['rhea_id'])
            reaction_smiles = row['reaction_smiles']
            ec_numbers = ec_mappings.get(rhea_id, [])

            reaction = self.reaction_smiles_to_reaction(
                rhea_id, reaction_smiles, chebi_names, ec_numbers
            )

            if reaction:
                reactions[rhea_id] = reaction
            else:
                failed_count += 1

        logger.info(f"Successfully loaded {len(reactions)} reactions from SMILES, {failed_count} failed")

        # Add polymer reactions from labels (no SMILES)
        if include_polymer_reactions and not rhea_ids:
            polymer_reactions = self.load_polymer_reactions()
            # Only add reactions with EC numbers for evaluation
            for rhea_id, reaction in polymer_reactions.items():
                if rhea_id in ec_mappings:
                    reactions[rhea_id] = reaction
            polymer_with_ec = sum(1 for rid in polymer_reactions if rid in ec_mappings)
            logger.info(f"Added {polymer_with_ec} polymer reactions with EC numbers")

        return reactions
        
    def apply_stoichiometry_from_label(
        self,
        reaction: Reaction,
        label: str,
        force: bool = False
    ) -> Reaction:
        """Apply stoichiometry and names from a parsed label to a reaction.

        Uses the label parser to extract (n) stoichiometry notation and
        applies it to matching participants in the reaction. Also fixes
        names when SMILES-based lookup gives wrong results (e.g., "S"
        mapping to hydrogen sulfide instead of sulfur atom).

        For polymer reactions with stoichiometry like "n X", the SMILES
        often has multiple copies (e.g., "S.S.S.S.S" for "5 sulfur").
        This method consolidates them into a single participant with
        the correct stoichiometry.

        Args:
            reaction: Reaction object with participants
            label: Human-readable reaction label with (n) notation
            force: If True, apply label matching even without (n) notation
                   (for polymer-type reactions like tRNA synthetases)

        Returns:
            Updated Reaction with stoichiometry/polymer_index/names fixed
        """
        if not force and not has_polymer_notation(label):
            return reaction

        left_parsed, right_parsed = parse_reaction_label(label)

        # Common molecule name aliases for better matching
        name_aliases = {
            'h2o': ['water', 'h2o'],
            'water': ['water', 'h2o'],
            'atp': ['atp', 'adenosine triphosphate'],
            'adp': ['adp', 'adenosine diphosphate'],
            'amp': ['amp', 'adenosine monophosphate'],
            'nad(+)': ['nad+', 'nad(+)', 'nicotinamide'],
            'nadh': ['nadh', 'reduced nad'],
            'nadp(+)': ['nadp+', 'nadp(+)'],
            'nadph': ['nadph'],
            'coa': ['coa', 'coenzyme a'],
            'h(+)': ['h+', 'h(+)', 'proton', 'hydrogen ion', 'hydron'],
            'h+': ['h+', 'h(+)', 'proton', 'hydron'],
            'hydron': ['h+', 'h(+)', 'proton', 'hydron'],
            'phosphate': ['phosphate', 'pi', 'orthophosphate'],
            'diphosphate': ['diphosphate', 'ppi', 'pyrophosphate'],
            'sulfur': ['sulfur', 's', 'compound_s', 'sulfur atom'],
            'hydrogen sulfide': ['hydrogen sulfide', 'h2s', 'hydrosulfide', 'sulfide'],
            'hydrosulfide': ['hydrogen sulfide', 'h2s', 'hydrosulfide', 'sulfide'],
            'h2': ['h2', 'hydrogen', 'molecular hydrogen', 'dihydrogen'],
            'dihydrogen': ['h2', 'hydrogen', 'molecular hydrogen', 'dihydrogen'],
            'o2': ['o2', 'oxygen', 'molecular oxygen', 'dioxygen'],
            'dioxygen': ['o2', 'oxygen', 'molecular oxygen', 'dioxygen'],
            'co2': ['co2', 'carbon dioxide'],
        }

        def names_match(
            name1: str, name2: str, parsed_polymer_type: Optional[str] = None
        ) -> bool:
            """Check if two names refer to the same molecule."""
            n1 = name1.lower().strip()
            n2 = name2.lower().strip()

            # Direct match
            if n1 == n2:
                return True

            # Partial match
            if n1 in n2 or n2 in n1:
                return True

            # Check aliases
            for key, aliases in name_aliases.items():
                if n1 == key or n1 in aliases:
                    if n2 == key or n2 in aliases:
                        return True

            # Handle compound_X names from SMILES
            if n1.startswith('compound_'):
                smiles_part = n1[9:]  # Remove "compound_"
                if smiles_part.lower() == n2 or n2 in smiles_part.lower():
                    return True

            # Handle "generic polymer" from SMILES matching to typed polymers from labels
            # If one side is "generic polymer" and the other has a polymer type, match them
            if 'generic polymer' in n1 or 'polymer' == n1:
                if parsed_polymer_type:
                    return True  # Any polymer type from label matches generic polymer

            return False

        def consolidate_and_match(
            participants: List[Participant], parsed_list
        ) -> List[Participant]:
            """Consolidate duplicate SMILES and match to parsed label participants.

            When SMILES has multiple copies (e.g., "S.S.S.S.S"), consolidate them
            into a single participant with the correct stoichiometry from the label.
            """
            if not parsed_list:
                return participants

            # Group participants by canonical SMILES
            from collections import defaultdict

            smiles_groups: Dict[str, List[Participant]] = defaultdict(list)
            non_smiles = []
            for p in participants:
                if p.smiles:
                    canonical = canonicalize_smiles(p.smiles)
                    if canonical:
                        smiles_groups[canonical].append(p)
                    else:
                        non_smiles.append(p)
                else:
                    non_smiles.append(p)

            # Build new participant list
            new_participants = []
            used_parsed = set()

            for canonical_smiles, group in smiles_groups.items():
                # If multiple copies, this might be a stoichiometric representation
                count = len(group)
                representative = group[0]

                # Try to match with a parsed participant
                matched_parsed = None
                for i, parsed in enumerate(parsed_list):
                    if i in used_parsed:
                        continue
                    # Match by name or by count pattern
                    if parsed.base_name:
                        # Check if parsed name matches (pass polymer_type for generic polymer matching)
                        if names_match(representative.name or "", parsed.base_name, parsed.polymer_type):
                            matched_parsed = parsed
                            used_parsed.add(i)
                            break
                        # If we have multiple copies and parsed has stoichiometry, likely a match
                        if count > 1 and parsed.stoichiometry:
                            # This could be the polymer species
                            matched_parsed = parsed
                            used_parsed.add(i)
                            break

                if matched_parsed:
                    # Use label name and stoichiometry
                    representative.name = matched_parsed.base_name
                    if matched_parsed.stoichiometry:
                        representative.stoichiometry = matched_parsed.stoichiometry
                    if matched_parsed.polymer_index:
                        representative.polymer_index = matched_parsed.polymer_index
                    if matched_parsed.polymer_type:
                        try:
                            representative.polymer_type = PolymerType(matched_parsed.polymer_type)
                        except ValueError:
                            representative.polymer_type = PolymerType.OTHER
                    if matched_parsed.linkage:
                        representative.linkage = matched_parsed.linkage
                    if matched_parsed.monomer:
                        representative.monomer = matched_parsed.monomer
                    if matched_parsed.location:
                        representative.location = matched_parsed.location
                    # Add only the representative (consolidated)
                    new_participants.append(representative)
                elif count > 1:
                    # Multiple copies but no match - consolidate anyway
                    representative.count = count
                    new_participants.append(representative)
                else:
                    # Single copy, keep all from group
                    new_participants.extend(group)

            # Handle non-SMILES participants - try to match "generic polymer" to parsed polymer types
            for p in non_smiles:
                if p.name and 'generic polymer' in p.name.lower():
                    # Try to find a matching parsed polymer by type
                    for i, parsed in enumerate(parsed_list):
                        if i in used_parsed:
                            continue
                        if parsed.polymer_type:
                            # Match generic polymer to typed parsed polymer
                            p.name = parsed.base_name
                            if parsed.polymer_type:
                                try:
                                    p.polymer_type = PolymerType(parsed.polymer_type)
                                except ValueError:
                                    p.polymer_type = PolymerType.OTHER
                            if parsed.stoichiometry:
                                p.stoichiometry = parsed.stoichiometry
                            if parsed.polymer_index:
                                p.polymer_index = parsed.polymer_index
                            used_parsed.add(i)
                            break
                new_participants.append(p)

            # Add any remaining parsed participants that weren't matched
            for i, parsed in enumerate(parsed_list):
                if i not in used_parsed and parsed.base_name:
                    # Create participant from label
                    new_p = Participant(
                        name=parsed.base_name,
                        stoichiometry=parsed.stoichiometry,
                        polymer_index=parsed.polymer_index,
                        location=parsed.location,
                    )
                    if parsed.polymer_type:
                        try:
                            new_p.polymer_type = PolymerType(parsed.polymer_type)
                        except ValueError:
                            new_p.polymer_type = PolymerType.OTHER
                    new_participants.append(new_p)

            return new_participants

        # Apply consolidation and matching
        reaction.left_participants = consolidate_and_match(
            reaction.left_participants, left_parsed
        )
        reaction.right_participants = consolidate_and_match(
            reaction.right_participants, right_parsed
        )

        return reaction

    def apply_locations_from_label(
        self,
        reaction: Reaction,
        label: str
    ) -> Reaction:
        """Apply location annotations from label to participants.

        For transport reactions, labels contain (in)/(out) annotations
        that indicate which side of a membrane a molecule is on.

        Args:
            reaction: Reaction object with participants
            label: Human-readable reaction label

        Returns:
            Updated Reaction with location annotations

        Examples:
            >>> etl = RheaTSVETL()
            >>> # Transport reaction: sulfate(out) + ATP + H2O = sulfate(in) + ADP + phosphate
        """
        # Only process if label has location annotations
        if '(in)' not in label and '(out)' not in label:
            return reaction

        left_parsed, right_parsed = parse_reaction_label(label)

        # Extract locations from parsed label by position
        # For transport reactions, the same molecule appears on both sides
        left_locations = {p.base_name: p.location for p in left_parsed if p.location}
        right_locations = {p.base_name: p.location for p in right_parsed if p.location}

        # Find ChEBI IDs shared between left and right (transported molecules)
        left_chebi = {p.chebi_id for p in reaction.left_participants if p.chebi_id}
        right_chebi = {p.chebi_id for p in reaction.right_participants if p.chebi_id}
        shared_chebi = left_chebi & right_chebi

        # For shared ChEBI IDs, assign locations based on which side they're on
        # Left side gets the location from left_parsed, right side from right_parsed
        if shared_chebi and (left_locations or right_locations):
            # Get any location from each side (transport has one annotated molecule per side)
            left_loc = next(iter(left_locations.values()), None)
            right_loc = next(iter(right_locations.values()), None)

            for participant in reaction.left_participants:
                if participant.chebi_id in shared_chebi and left_loc:
                    participant.location = left_loc

            for participant in reaction.right_participants:
                if participant.chebi_id in shared_chebi and right_loc:
                    participant.location = right_loc

        return reaction

    def create_rhea_terms(
        self,
        reactions: Dict[str, Reaction],
        include_ec: bool = True,
        include_labels: bool = True,
        include_go: bool = True,
        go_cache_path: Optional[Path] = None,
        master_only: bool = True
    ) -> List[RheaTerm]:
        """Convert reactions to RheaTerm objects for evaluation.

        Args:
            reactions: Dictionary of RHEA ID to Reaction
            include_ec: Whether to include EC number information
            include_labels: Whether to load and apply human-readable labels
            include_go: Whether to include GO term mappings
            go_cache_path: Path to go_terms.jsonl file
            master_only: If True, consolidate to master IDs only (direction-neutral)

        Returns:
            List of RheaTerm objects
        """
        ec_mappings = self.load_ec_mappings() if include_ec else {}
        go_mappings = self.load_go_mappings(go_cache_path) if include_go else {}

        # Load labels for all reactions (includes polymer notation)
        labels = {}
        if include_labels:
            labels = self.load_reaction_labels()

        # Load EQUATION ChEBI lookup for normalization
        equation_lookup = self.load_equation_chebi_lookup()

        # Load directions file for proper master ID mapping
        # (the arithmetic formula (id // 4) * 4 is unreliable)
        directions = self.load_directions()

        rhea_terms = []
        polymer_count = 0
        chebi_normalized_count = 0
        seen_master_ids = set()

        for rhea_id, reaction in reactions.items():
            # Convert to master ID using directions file (more reliable than arithmetic)
            master_id = directions.get(rhea_id, rhea_id)

            if master_only:
                if master_id in seen_master_ids:
                    continue  # Skip duplicate master IDs
                seen_master_ids.add(master_id)
                display_id = master_id
            else:
                display_id = rhea_id

            # Get label - try master ID first, then original
            label = labels.get(master_id) or labels.get(rhea_id, f"Reaction {display_id}")

            # Apply stoichiometry/polymer info from label if needed
            # Check for (n) notation OR if reaction has "generic polymer" participants
            has_generic_polymer = any(
                (p.name and 'generic polymer' in p.name.lower())
                for p in reaction.left_participants + reaction.right_participants
            )
            if has_polymer_notation(label) or has_generic_polymer:
                # Force label matching for generic polymer reactions even without (n) notation
                reaction = self.apply_stoichiometry_from_label(reaction, label, force=has_generic_polymer)
                polymer_count += 1

            # Normalize ChEBI IDs using EQUATION data and name-based lookup
            # MUST happen BEFORE location application (locations need shared ChEBI)
            # Count participants without ChEBI before normalization
            missing_before = sum(
                1 for p in reaction.left_participants + reaction.right_participants
                if not p.chebi_id
            )

            # Apply normalization using the ORIGINAL rhea_id for EQUATION lookup
            # (not display_id/master_id, which may map to a different reaction)
            self.normalize_reaction_chebi(reaction, rhea_id, equation_lookup)

            # Count how many were filled in
            missing_after = sum(
                1 for p in reaction.left_participants + reaction.right_participants
                if not p.chebi_id
            )
            if missing_before > missing_after:
                chebi_normalized_count += (missing_before - missing_after)

            # Apply location annotations from label (for transport reactions)
            # Must happen AFTER ChEBI normalization so shared_chebi detection works
            self.apply_locations_from_label(reaction, label)

            # Get GO terms for master ID (go_mappings uses int keys)
            go_terms = go_mappings.get(int(master_id), [])

            # Get EC numbers - try master ID first
            ec_numbers = ec_mappings.get(master_id, []) or ec_mappings.get(rhea_id, [])

            # Create RheaTerm
            rhea_term = RheaTerm(
                rhea_id=f"RHEA:{display_id}",
                label=label,
                reaction=reaction,
                ec_numbers=ec_numbers,
                go_terms=go_terms,
                direction="bidirectional",
                parent_rhea_id=None
            )
            rhea_terms.append(rhea_term)

        logger.info(
            f"Created {len(rhea_terms)} RheaTerm objects "
            f"({polymer_count} with polymer notation, "
            f"{chebi_normalized_count} ChEBI IDs normalized)"
        )
        return rhea_terms


# Convenience functions for easy usage
def load_rhea_reactions_tsv(
    limit: Optional[int] = None,
    cache_dir: Optional[Path] = None,
    go_cache_path: Optional[Path] = None
) -> Dict[str, Reaction]:
    """Load RHEA reactions using TSV approach.

    Args:
        limit: Optional limit on reactions to load
        cache_dir: Optional cache directory
        go_cache_path: Path to GO cache (for including GO-mapped reactions)

    Returns:
        Dictionary of RHEA ID to Reaction objects
    """
    if go_cache_path is None:
        go_cache_path = Path("cache") / "go_terms.jsonl"
    etl = RheaTSVETL(cache_dir=cache_dir or Path("cache/rhea_tsv"))
    return etl.load_reactions_batch(limit=limit, go_cache_path=go_cache_path)


def load_rhea_terms_tsv(
    limit: Optional[int] = None,
    cache_dir: Optional[Path] = None,
    go_cache_path: Optional[Path] = None
) -> List[RheaTerm]:
    """Load RHEA terms using TSV approach for evaluation.

    Args:
        limit: Optional limit on reactions to load
        cache_dir: Optional cache directory
        go_cache_path: Path to GO cache (for including GO-mapped reactions)

    Returns:
        List of RheaTerm objects
    """
    if go_cache_path is None:
        go_cache_path = Path("cache") / "go_terms.jsonl"
    etl = RheaTSVETL(cache_dir=cache_dir or Path("cache/rhea_tsv"))
    reactions = etl.load_reactions_batch(limit=limit, go_cache_path=go_cache_path)
    return etl.create_rhea_terms(reactions, go_cache_path=go_cache_path)


if __name__ == "__main__":
    # Demo usage
    logging.basicConfig(level=logging.INFO)
    
    print("=== RHEA TSV ETL Demo ===")
    
    # Load a small sample
    reactions = load_rhea_reactions_tsv(limit=10)
    
    print(f"Loaded {len(reactions)} reactions")
    
    # Show first reaction
    if reactions:
        first_id = next(iter(reactions.keys()))
        first_reaction = reactions[first_id]
        print(f"\nSample reaction {first_id}:")
        print(f"  Left: {len(first_reaction.left_participants)} participants")
        print(f"  Right: {len(first_reaction.right_participants)} participants")
        
        for i, p in enumerate(first_reaction.left_participants[:3]):
            print(f"    L{i}: {p.name} ({p.smiles})")
