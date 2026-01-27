"""
Enhanced Redock Analysis Script

Properly handles actives (self-docking) vs decoys (enrichment) analysis.

Usage:
    python redock_analysis_fixed.py results.csv --output report.json
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse
from dataclasses import dataclass, asdict
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


@dataclass
class ActiveStatistics:
    """Statistics for active (self-docking) cases."""
    n_cases: int
    n_valid_rmsd: int
    success_rate: float
    mean_rmsd: Optional[float]
    median_rmsd: Optional[float]
    min_rmsd: Optional[float]
    max_rmsd: Optional[float]
    mean_runtime: float
    

@dataclass
class EnrichmentStatistics:
    """Statistics for enrichment (active vs decoy discrimination)."""
    n_actives: int
    n_decoys: int
    roc_auc: Optional[float]
    ef_1_percent: Optional[float]
    ef_5_percent: Optional[float]
    ef_10_percent: Optional[float]
    mean_active_score: float
    mean_decoy_score: float
    score_separation: float


class RedockAnalyzer:
    """
    Enhanced redock analyzer that properly separates actives and decoys.
    """
    
    def __init__(
        self,
        rmsd_threshold: float = 2.0,
        active_similarity_threshold: float = 0.85
    ):
        """
        Initialize analyzer.
        
        Args:
            rmsd_threshold: RMSD threshold for success (Angstroms)
            active_similarity_threshold: Similarity threshold for active detection
        """
        self.rmsd_threshold = rmsd_threshold
        self.active_similarity_threshold = active_similarity_threshold
    
    def classify_results(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Classify each result as active or decoy based on RMSD availability.
        
        Args:
            df: DataFrame with redock results
            
        Returns:
            DataFrame with 'molecule_type' column added
        """
        df = df.copy()
        
        # Strategy: If RMSD = 999.9, it's likely a decoy (molecule mismatch)
        # If RMSD < 999, it's likely an active (same molecule)
        def classify_row(row):
            rmsd = row.get('best_rmsd', 999.9)
            
            # Check if RMSD calculation failed (999.9 is sentinel)
            if pd.isna(rmsd) or rmsd >= 999:
                return 'decoy'
            else:
                return 'active'
        
        df['molecule_type'] = df.apply(classify_row, axis=1)
        
        # Override with explicit control labels when present
        if 'control_label' in df.columns:
            # numeric labels
            df.loc[df['control_label'] == 0, 'molecule_type'] = 'decoy'
            df.loc[df['control_label'] == 1, 'molecule_type'] = 'active'
            # string labels
            df.loc[df['control_label'] == 'decoy', 'molecule_type'] = 'decoy'
            df.loc[df['control_label'] == 'active', 'molecule_type'] = 'active'
        
        return df
    
    def analyze_actives(self, df: pd.DataFrame) -> ActiveStatistics:
        """
        Analyze active (self-docking) cases only.
        
        Args:
            df: DataFrame with redock results
            
        Returns:
            ActiveStatistics object
        """
        # Filter to actives with valid RMSD
        actives = df[df['molecule_type'] == 'active'].copy()
        valid_actives = actives[actives['best_rmsd'] < 999].copy()
        
        if len(valid_actives) == 0:
            logger.warning("No valid active cases found")
            return ActiveStatistics(
                n_cases=len(actives),
                n_valid_rmsd=0,
                success_rate=0.0,
                mean_rmsd=None,
                median_rmsd=None,
                min_rmsd=None,
                max_rmsd=None,
                mean_runtime=0.0
            )
        
        # Calculate success rate
        valid_actives['success'] = valid_actives['best_rmsd'] < self.rmsd_threshold
        success_rate = valid_actives['success'].mean() * 100
        
        return ActiveStatistics(
            n_cases=len(actives),
            n_valid_rmsd=len(valid_actives),
            success_rate=success_rate,
            mean_rmsd=valid_actives['best_rmsd'].mean(),
            median_rmsd=valid_actives['best_rmsd'].median(),
            min_rmsd=valid_actives['best_rmsd'].min(),
            max_rmsd=valid_actives['best_rmsd'].max(),
            mean_runtime=valid_actives['runtime_sec'].mean()
        )
    
    def analyze_enrichment(self, df: pd.DataFrame) -> EnrichmentStatistics:
        """
        Analyze enrichment (active vs decoy discrimination).
        
        Args:
            df: DataFrame with redock results
            
        Returns:
            EnrichmentStatistics object
        """
        # Only include explicit controls (labelled actives/decoys)
        if 'control_label' in df.columns:
            control_df = df[df['control_label'].isin([0, 1])].copy()
            actives = control_df[control_df['control_label'] == 1].copy()
            decoys = control_df[control_df['control_label'] == 0].copy()
        else:
            actives = df[df['molecule_type'] == 'active'].copy()
            decoys = df[df['molecule_type'] == 'decoy'].copy()
        
        if len(actives) == 0 or len(decoys) == 0:
            logger.warning("Need both actives and decoys for enrichment analysis")
            return EnrichmentStatistics(
                n_actives=len(actives),
                n_decoys=len(decoys),
                roc_auc=None,
                ef_1_percent=None,
                ef_5_percent=None,
                ef_10_percent=None,
                mean_active_score=0.0,
                mean_decoy_score=0.0,
                score_separation=0.0
            )
        
        # Prepare data for ROC
        labels = [1] * len(actives) + [0] * len(decoys)
        
        # Use best score (lower is better, so negate for ROC)
        active_scores = actives['best_score'].values
        decoy_scores = decoys['best_score'].values
        scores = list(-active_scores) + list(-decoy_scores)  # Negate for ROC
        
        # Calculate ROC AUC
        try:
            from sklearn.metrics import roc_auc_score
            roc_auc = roc_auc_score(labels, scores)
        except Exception as e:
            logger.warning(f"ROC AUC calculation failed: {e}")
            roc_auc = None
        
        # Calculate enrichment factors
        def calculate_ef(percent: float) -> Optional[float]:
            try:
                n_total = len(labels)
                n_actives = sum(labels)
                n_select = int(n_total * percent / 100)
                
                if n_select == 0:
                    return None
                
                # Sort by score (descending = better)
                sorted_indices = sorted(
                    range(len(scores)),
                    key=lambda i: scores[i],
                    reverse=True
                )
                
                # Count actives in top N
                actives_found = sum(labels[i] for i in sorted_indices[:n_select])
                
                # EF = (actives_found / n_select) / (n_actives / n_total)
                if n_select == 0:
                    return None
                ef = (actives_found / n_select) / (n_actives / n_total)
                return ef
                
            except Exception as e:
                logger.warning(f"EF calculation failed: {e}")
                return None
        
        # Score statistics
        mean_active_score = active_scores.mean()
        mean_decoy_score = decoy_scores.mean()
        score_separation = mean_decoy_score - mean_active_score  # Should be positive
        
        return EnrichmentStatistics(
            n_actives=len(actives),
            n_decoys=len(decoys),
            roc_auc=roc_auc,
            ef_1_percent=calculate_ef(1.0),
            ef_5_percent=calculate_ef(5.0),
            ef_10_percent=calculate_ef(10.0),
            mean_active_score=mean_active_score,
            mean_decoy_score=mean_decoy_score,
            score_separation=score_separation
        )
    
    def generate_report(self, df: pd.DataFrame) -> Dict:
        """
        Generate comprehensive report separating active and enrichment analysis.
        
        Args:
            df: DataFrame with redock results
            
        Returns:
            Dictionary with complete analysis
        """
        # Classify results
        df = self.classify_results(df)
        
        # Analyze actives
        active_stats = self.analyze_actives(df)
        
        # Analyze enrichment
        enrichment_stats = self.analyze_enrichment(df)
        
        # Overall statistics
        report = {
            'summary': {
                'total_cases': len(df),
                'n_actives': len(df[df['molecule_type'] == 'active']),
                'n_decoys': len(df[df['molecule_type'] == 'decoy']),
                'rmsd_threshold': self.rmsd_threshold
            },
            'active_analysis': asdict(active_stats),
            'enrichment_analysis': asdict(enrichment_stats),
            'interpretation': self._generate_interpretation(
                active_stats,
                enrichment_stats
            )
        }
        
        return report
    
    def _generate_interpretation(
        self,
        active_stats: ActiveStatistics,
        enrichment_stats: EnrichmentStatistics
    ) -> Dict:
        """Generate human-readable interpretation of results."""
        interpretation = {}
        
        # Active performance
        if active_stats.n_valid_rmsd > 0:
            if active_stats.success_rate >= 70:
                interpretation['active_quality'] = 'Excellent'
                interpretation['active_message'] = (
                    f"Self-docking performs excellently with "
                    f"{active_stats.success_rate:.1f}% success rate"
                )
            elif active_stats.success_rate >= 50:
                interpretation['active_quality'] = 'Good'
                interpretation['active_message'] = (
                    f"Self-docking performs well with "
                    f"{active_stats.success_rate:.1f}% success rate"
                )
            elif active_stats.success_rate >= 30:
                interpretation['active_quality'] = 'Fair'
                interpretation['active_message'] = (
                    f"Self-docking shows moderate performance with "
                    f"{active_stats.success_rate:.1f}% success rate"
                )
            else:
                interpretation['active_quality'] = 'Poor'
                interpretation['active_message'] = (
                    f"Self-docking needs improvement with only "
                    f"{active_stats.success_rate:.1f}% success rate"
                )
        else:
            interpretation['active_quality'] = 'Unknown'
            interpretation['active_message'] = 'No valid self-docking cases'
        
        # Enrichment performance
        if enrichment_stats.roc_auc is not None:
            if enrichment_stats.roc_auc >= 0.7:
                interpretation['enrichment_quality'] = 'Excellent'
                interpretation['enrichment_message'] = (
                    f"Excellent discrimination between actives and decoys "
                    f"(ROC AUC = {enrichment_stats.roc_auc:.3f})"
                )
            elif enrichment_stats.roc_auc >= 0.6:
                interpretation['enrichment_quality'] = 'Good'
                interpretation['enrichment_message'] = (
                    f"Good discrimination between actives and decoys "
                    f"(ROC AUC = {enrichment_stats.roc_auc:.3f})"
                )
            elif enrichment_stats.roc_auc >= 0.5:
                interpretation['enrichment_quality'] = 'Fair'
                interpretation['enrichment_message'] = (
                    f"Modest discrimination between actives and decoys "
                    f"(ROC AUC = {enrichment_stats.roc_auc:.3f})"
                )
            else:
                interpretation['enrichment_quality'] = 'Poor'
                interpretation['enrichment_message'] = (
                    f"Poor discrimination between actives and decoys "
                    f"(ROC AUC = {enrichment_stats.roc_auc:.3f}). "
                    f"May need parameter optimization."
                )
        else:
            interpretation['enrichment_quality'] = 'Unknown'
            interpretation['enrichment_message'] = 'Insufficient data for enrichment analysis'
        
        return interpretation
    
    def export_report(
        self,
        df: pd.DataFrame,
        output_json: Optional[Path] = None,
        output_markdown: Optional[Path] = None
    ) -> Dict:
        """
        Generate and export report.
        
        Args:
            df: DataFrame with results
            output_json: Path for JSON report (optional)
            output_markdown: Path for Markdown report (optional)
            
        Returns:
            Report dictionary
        """
        report = self.generate_report(df)
        
        # Export JSON
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"JSON report saved: {output_json}")
        
        # Export Markdown
        if output_markdown:
            md = self._format_markdown(report)
            with open(output_markdown, 'w') as f:
                f.write(md)
            logger.info(f"Markdown report saved: {output_markdown}")
        
        return report
    
    def _format_markdown(self, report: Dict) -> str:
        """Format report as Markdown."""
        lines = []
        
        lines.append("# Redock Analysis Report (Enhanced)")
        lines.append("")
        
        # Summary
        summary = report['summary']
        lines.append("## Summary")
        lines.append(f"- Total cases: {summary['total_cases']}")
        lines.append(f"- Active cases: {summary['n_actives']}")
        lines.append(f"- Decoy cases: {summary['n_decoys']}")
        lines.append(f"- RMSD threshold: {summary['rmsd_threshold']} Å")
        lines.append("")
        
        # Active Analysis
        active = report['active_analysis']
        lines.append("## Self-Docking Performance (Actives Only)")
        lines.append(f"- Valid cases: {active['n_valid_rmsd']} / {active['n_cases']}")
        lines.append(f"- **Success rate: {active['success_rate']:.1f}%**")
        if active['mean_rmsd'] is not None:
            lines.append(f"- Mean RMSD: {active['mean_rmsd']:.2f} Å")
            lines.append(f"- Median RMSD: {active['median_rmsd']:.2f} Å")
            lines.append(f"- RMSD range: [{active['min_rmsd']:.2f}, {active['max_rmsd']:.2f}] Å")
        lines.append(f"- Mean runtime: {active['mean_runtime']:.2f} sec")
        lines.append("")
        
        # Enrichment Analysis
        enrich = report['enrichment_analysis']
        lines.append("## Enrichment Performance (Active vs Decoy)")
        lines.append(f"- Actives: {enrich['n_actives']}")
        lines.append(f"- Decoys: {enrich['n_decoys']}")
        if enrich['roc_auc'] is not None:
            lines.append(f"- **ROC AUC: {enrich['roc_auc']:.3f}**")
        if enrich['ef_1_percent'] is not None:
            lines.append(f"- EF @ 1%: {enrich['ef_1_percent']:.2f}x")
        if enrich['ef_5_percent'] is not None:
            lines.append(f"- EF @ 5%: {enrich['ef_5_percent']:.2f}x")
        if enrich['ef_10_percent'] is not None:
            lines.append(f"- EF @ 10%: {enrich['ef_10_percent']:.2f}x")
        lines.append(f"- Mean active score: {enrich['mean_active_score']:.2f}")
        lines.append(f"- Mean decoy score: {enrich['mean_decoy_score']:.2f}")
        lines.append(f"- Score separation: {enrich['score_separation']:.2f}")
        lines.append("")
        
        # Interpretation
        interp = report['interpretation']
        lines.append("## Interpretation")
        lines.append(f"- **Self-docking: {interp.get('active_quality', 'Unknown')}**")
        lines.append(f"  - {interp.get('active_message', 'N/A')}")
        lines.append(f"- **Enrichment: {interp.get('enrichment_quality', 'Unknown')}**")
        lines.append(f"  - {interp.get('enrichment_message', 'N/A')}")
        lines.append("")
        
        return "\n".join(lines)


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Enhanced redock analysis with proper active/decoy handling"
    )
    parser.add_argument(
        'input_csv',
        type=Path,
        help='Input CSV file with redock results'
    )
    parser.add_argument(
        '--output-json',
        type=Path,
        default=None,
        help='Output path for JSON report'
    )
    parser.add_argument(
        '--output-md',
        type=Path,
        default=None,
        help='Output path for Markdown report'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=2.0,
        help='RMSD threshold for success (default: 2.0 Å)'
    )
    
    args = parser.parse_args()
    
    # Load results
    logger.info(f"Loading results from {args.input_csv}")
    df = pd.read_csv(args.input_csv)
    
    # Run analysis
    analyzer = RedockAnalyzer(rmsd_threshold=args.threshold)
    report = analyzer.export_report(
        df,
        output_json=args.output_json,
        output_markdown=args.output_md
    )
    
    # Print summary
    print("\n" + "="*60)
    print("REDOCK ANALYSIS SUMMARY")
    print("="*60)
    
    active = report['active_analysis']
    print(f"\nSelf-Docking (Actives):")
    print(f"  Success Rate: {active['success_rate']:.1f}%")
    if active['mean_rmsd'] is not None:
        print(f"  Mean RMSD: {active['mean_rmsd']:.2f} Å")
    
    enrich = report['enrichment_analysis']
    print(f"\nEnrichment (Active vs Decoy):")
    if enrich['roc_auc'] is not None:
        print(f"  ROC AUC: {enrich['roc_auc']:.3f}")
    print(f"  Actives: {enrich['n_actives']}, Decoys: {enrich['n_decoys']}")
    
    interp = report['interpretation']
    print(f"\nInterpretation:")
    print(f"  {interp.get('active_message', 'N/A')}")
    print(f"  {interp.get('enrichment_message', 'N/A')}")
    print()


if __name__ == '__main__':
    main()
