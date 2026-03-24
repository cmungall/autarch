"""Comprehensive summary of reaction classes with GO, EC, and RHEA information."""

from typing import List, Optional, Type
from dataclasses import dataclass

from oaklib import get_adapter  # type: ignore[import-untyped]
from rich.console import Console
from rich.table import Table

from autarch.ontology.reaction import ReactionClass


@dataclass
class ClassSummary:
    """Summary information for a reaction class."""
    class_name: str
    go_id: Optional[str]
    go_label: Optional[str]
    ec_prefix: Optional[str]
    ec_mappings: List[str]
    rhea_mappings: List[str]


class ClassSummarizer:
    """Generates comprehensive summaries of all reaction classes."""
    
    def __init__(self):
        """Initialize with GO ontology adapter."""
        self.go_adapter = get_adapter("sqlite:obo:go")
        self.console = Console()
        
    def get_go_info(self, go_id: str) -> Optional[str]:
        """Get GO term label."""
        try:
            return self.go_adapter.label(go_id)
        except Exception:
            return None
            
    def get_ec_mappings(self, go_id: str) -> List[str]:
        """Get EC mappings for a GO term - simplified version for speed."""
        # For now, return empty to keep it fast
        # TODO: Implement efficient EC mapping retrieval
        return []
            
    def get_rhea_mappings(self, go_id: str) -> List[str]:
        """Get RHEA mappings for a GO term - simplified version for speed."""
        # For now, return empty to keep it fast  
        # TODO: Implement efficient RHEA mapping retrieval
        return []
            
    def summarize_class(self, cls: Type[ReactionClass]) -> ClassSummary:
        """Create summary for a single reaction class."""
        class_name = cls.__name__
        go_id = getattr(cls, 'GO_ID', None)
        ec_prefix = getattr(cls, 'EC_NUMBER_PREFIX', None)
        
        go_label = None
        ec_mappings = []
        rhea_mappings = []
        
        if go_id:
            go_label = self.get_go_info(go_id)
            ec_mappings = self.get_ec_mappings(go_id)
            rhea_mappings = self.get_rhea_mappings(go_id)
            
        return ClassSummary(
            class_name=class_name,
            go_id=go_id,
            go_label=go_label,
            ec_prefix=ec_prefix,
            ec_mappings=ec_mappings,
            rhea_mappings=rhea_mappings
        )
        
    def summarize_all_classes(self) -> List[ClassSummary]:
        """Summarize all reaction classes."""
        summaries = []
        
        # Get all reaction classes from ontology module
        import autarch.ontology as ont_module
        
        for name in dir(ont_module):
            cls = getattr(ont_module, name)
            
            # Check if it's a ReactionClass
            if (hasattr(cls, 'GO_ID') and 
                hasattr(cls, '__bases__') and 
                isinstance(cls, type)):
                
                summary = self.summarize_class(cls)
                summaries.append(summary)
                
        # Sort by class name
        return sorted(summaries, key=lambda x: x.class_name)
        
    def print_summary_table(self, summaries: List[ClassSummary]):
        """Print a rich table of class summaries."""
        table = Table(title="Autarch Reaction Class Summary")
        
        table.add_column("Class Name", style="cyan", width=15, no_wrap=True)
        table.add_column("GO ID", style="green", width=12, no_wrap=True)  
        table.add_column("GO Label", style="green", width=25)
        table.add_column("EC Prefix", style="yellow", width=9, no_wrap=True)
        table.add_column("RHEA", style="blue", width=5, no_wrap=True)
        table.add_column("Notes", style="magenta", width=18)
        
        for summary in summaries:
            go_id = summary.go_id or "None"
            go_label = summary.go_label or "No label found"
            ec_prefix = summary.ec_prefix or "None"
            
            # Add notes for special classes
            notes = []
            if "Protein" in summary.class_name and "Ligase" in summary.class_name:
                notes.append("Ultra-specific")
            if summary.class_name in ["EpoxideHydrolase", "CytochromeCOxidase", "Methyltransferase"]:
                notes.append("Perfect performance")
            if summary.class_name.startswith("Carbon") or summary.class_name.startswith("SUMO"):
                notes.append("GO-compliant naming") 
                
            notes_str = ", ".join(notes) if notes else ""
            
            # Show full GO ID
            display_go_id = go_id if go_id != "None" else "-"
            
            # Truncate long labels
            display_label = go_label[:25] + "..." if len(go_label) > 28 else go_label
            
            # RHEA count placeholder (TODO: implement)
            rhea_count = str(len(summary.rhea_mappings)) if summary.rhea_mappings else "-"
            
            table.add_row(
                summary.class_name,
                display_go_id,
                display_label,
                ec_prefix or "-",
                rhea_count,
                notes_str
            )
            
        self.console.print(table)
        
        # Print summary statistics
        total_classes = len(summaries)
        with_go = sum(1 for s in summaries if s.go_id)
        with_ec_prefix = sum(1 for s in summaries if s.ec_prefix)
        ultra_specific = sum(1 for s in summaries if "Protein" in s.class_name and "Ligase" in s.class_name)
        new_classes = [
            "BiotinBiotinCarboxylCarrierProteinLigase",
            "LigaseFormingCarbonNitrogenBonds",
            "CytochromeCOxidase",
            "EpoxideHydrolase",
            "UbiquitinProteinLigase",
            "Methyltransferase",
            "NADHOrNADPHDehydrogenaseQuinone",
            "SUMOTransferase",
            "SqualeneMonooxygenase",
        ]
        new_count = sum(1 for s in summaries if s.class_name in new_classes)
        
        self.console.print("\n📊 Summary Statistics:")
        self.console.print(f"  • Total Classes: {total_classes}")
        self.console.print(f"  • With GO ID: {with_go}")
        self.console.print(f"  • With EC Prefix: {with_ec_prefix}")
        self.console.print(f"  • Ultra-Specific Classes: {ultra_specific}")
        self.console.print(f"  • Newly Added Classes: {new_count}")
        
        # Highlight some key classes
        self.console.print("\n🌟 Key New Ultra-Specific Classes:")
        for s in summaries:
            if s.class_name in [
                "BiotinBiotinCarboxylCarrierProteinLigase",
                "EpoxideHydrolase",
                "CytochromeCOxidase",
                "SUMOTransferase",
            ]:
                self.console.print(f"  • {s.class_name}: {s.go_label}")
        
    def generate_summary(self):
        """Generate and print the complete summary."""
        self.console.print("🔍 Generating comprehensive class summary...\n")
        
        summaries = self.summarize_all_classes()
        self.print_summary_table(summaries)
