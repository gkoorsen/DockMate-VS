"""
Ensemble validation experiment runner.
Compares single structure vs ensemble docking performance.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from datetime import datetime
from loguru import logger
import matplotlib.pyplot as plt
import seaborn as sns

from ..analysis.enrichment_analysis import EnrichmentAnalysis, BootstrapValidation


class EnsembleValidationExperiment:
    """
    Run comprehensive validation of ensemble vs single structure docking.
    
    Workflow:
    1. Dock all compounds to single best structure
    2. Dock all compounds to ensemble (5 structures)
    3. Calculate enrichment metrics for both
    4. Statistical comparison with bootstrap
    5. Generate visualizations and report
    """
    
    def __init__(self, docking_manager, output_dir='ensemble_validation'):
        """
        Initialize experiment.
        
        Parameters
        ----------
        docking_manager : object
            Docking manager with methods:
            - dock_compound(compound, protein, grid_params) -> dict
            - dock_ensemble(compound, proteins, grid_params) -> dict
        output_dir : str or Path
            Directory for outputs
        """
        self.docking_manager = docking_manager
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Results storage
        self.single_results = None
        self.ensemble_results = None
        self.comparison = None
        
        logger.info(f"Experiment initialized. Output: {self.output_dir}")
    
    def load_benchmark_dataset(self, actives_file, decoys_file):
        """
        Load actives and decoys for benchmarking.
        
        Parameters
        ----------
        actives_file : str or Path
            File with active compounds (SDF, CSV, etc.)
        decoys_file : str or Path
            File with decoy compounds
        
        Returns
        -------
        pd.DataFrame : Combined dataset with columns:
            - compound_id
            - compound_data (mol object or SMILES)
            - is_active (bool)
        """
        logger.info("Loading benchmark dataset...")
        
        # Load actives
        actives_df = self._load_compound_file(actives_file, is_active=True)
        logger.info(f"  Loaded {len(actives_df)} actives from {actives_file}")
        
        # Load decoys
        decoys_df = self._load_compound_file(decoys_file, is_active=False)
        logger.info(f"  Loaded {len(decoys_df)} decoys from {decoys_file}")
        
        # Combine
        dataset = pd.concat([actives_df, decoys_df], ignore_index=True)
        
        ratio = len(decoys_df) / len(actives_df)
        logger.info(f"  Total: {len(dataset)} compounds (ratio 1:{ratio:.1f})")
        
        return dataset
    
    def _load_compound_file(self, filepath, is_active):
        """Load compounds from file (SDF, CSV, etc.)."""
        filepath = Path(filepath)
        
        if filepath.suffix == '.sdf':
            return self._load_sdf(filepath, is_active)
        elif filepath.suffix == '.csv':
            return self._load_csv(filepath, is_active)
        else:
            raise ValueError(f"Unsupported file format: {filepath.suffix}")
    
    def _load_sdf(self, filepath, is_active):
        """Load SDF file."""
        from rdkit import Chem
        
        compounds = []
        supplier = Chem.SDMolSupplier(str(filepath))
        
        for i, mol in enumerate(supplier):
            if mol is None:
                logger.warning(f"Could not read molecule {i} from {filepath}")
                continue
            
            # Get compound ID
            if mol.HasProp('_Name'):
                compound_id = mol.GetProp('_Name')
            else:
                compound_id = f"compound_{i}"
            
            compounds.append({
                'compound_id': compound_id,
                'compound_data': mol,
                'is_active': is_active
            })
        
        return pd.DataFrame(compounds)
    
    def _load_csv(self, filepath, is_active):
        """Load CSV file (expects columns: compound_id, smiles)."""
        from rdkit import Chem
        
        df = pd.read_csv(filepath)
        
        if 'smiles' not in df.columns or 'compound_id' not in df.columns:
            raise ValueError("CSV must have 'compound_id' and 'smiles' columns")
        
        # Convert SMILES to mol objects
        df['compound_data'] = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x))
        df['is_active'] = is_active
        
        # Remove invalid molecules
        df = df[df['compound_data'].notna()].reset_index(drop=True)
        
        return df[['compound_id', 'compound_data', 'is_active']]
    
    def run_full_experiment(
        self,
        benchmark_dataset,
        single_structure,
        ensemble_structures,
        grid_params,
        n_bootstrap=1000
    ):
        """
        Run complete validation experiment.
        
        Parameters
        ----------
        benchmark_dataset : pd.DataFrame
            Dataset from load_benchmark_dataset()
        single_structure : str or Path
            Path to single best structure
        ensemble_structures : list of str/Path
            Paths to all ensemble structures
        grid_params : dict
            Grid box parameters (center, size)
        n_bootstrap : int
            Number of bootstrap iterations
        
        Returns
        -------
        dict : Complete results
        """
        start_time = datetime.now()
        
        logger.info("\n" + "="*70)
        logger.info("ENSEMBLE VALIDATION EXPERIMENT")
        logger.info("="*70)
        logger.info(f"Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"Dataset: {len(benchmark_dataset)} compounds")
        logger.info(f"  Actives: {benchmark_dataset['is_active'].sum()}")
        logger.info(f"  Decoys: {(~benchmark_dataset['is_active']).sum()}")
        logger.info(f"Single structure: {single_structure}")
        logger.info(f"Ensemble: {len(ensemble_structures)} structures")
        
        # Run experiments
        results = {}
        
        # 1. Single structure docking
        logger.info("\n" + "-"*70)
        logger.info("[1/5] SINGLE STRUCTURE DOCKING")
        logger.info("-"*70)
        self.single_results = self._run_single_experiment(
            benchmark_dataset,
            single_structure,
            grid_params
        )
        results['single_results'] = self.single_results
        
        # 2. Ensemble docking
        logger.info("\n" + "-"*70)
        logger.info("[2/5] ENSEMBLE DOCKING")
        logger.info("-"*70)
        self.ensemble_results = self._run_ensemble_experiment(
            benchmark_dataset,
            ensemble_structures,
            grid_params
        )
        results['ensemble_results'] = self.ensemble_results
        
        # 3. Calculate enrichment metrics
        logger.info("\n" + "-"*70)
        logger.info("[3/5] ENRICHMENT ANALYSIS")
        logger.info("-"*70)
        
        single_analysis = EnrichmentAnalysis(self.single_results)
        ensemble_analysis = EnrichmentAnalysis(self.ensemble_results)
        
        single_metrics = single_analysis.calculate_all_metrics()
        ensemble_metrics = ensemble_analysis.calculate_all_metrics()
        
        results['single_metrics'] = single_metrics
        results['ensemble_metrics'] = ensemble_metrics
        
        # Print summaries
        single_analysis.print_summary("Single Structure")
        ensemble_analysis.print_summary("Ensemble")
        
        # 4. Statistical comparison
        logger.info("-"*70)
        logger.info("[4/5] STATISTICAL VALIDATION")
        logger.info("-"*70)
        
        bootstrap_results = BootstrapValidation.compare_methods(
            self.single_results,
            self.ensemble_results,
            n_iterations=n_bootstrap
        )
        results['bootstrap'] = bootstrap_results
        
        BootstrapValidation.print_bootstrap_results(
            bootstrap_results,
            method_names=('Single', 'Ensemble')
        )

        # Make recommendation and save core results before plotting.
        logger.info("\n" + "="*70)
        logger.info("RECOMMENDATION")
        logger.info("="*70)
        recommendation = self._make_recommendation(results)
        results['recommendation'] = recommendation
        self._save_results(results)

        # 5. Generate visualizations (non-blocking for saved results)
        logger.info("-"*70)
        logger.info("[5/5] GENERATING VISUALIZATIONS")
        logger.info("-"*70)
        try:
            self._generate_all_plots(
                single_analysis,
                ensemble_analysis,
                bootstrap_results
            )
        except Exception as exc:
            logger.warning("Plot generation failed: {}", exc)
        
        # Timing
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        logger.info(f"\nCompleted in {duration:.1f} seconds")
        logger.info(f"Results saved to: {self.output_dir}")
        
        return results
    
    def _run_single_experiment(self, dataset, structure, grid_params):
        """Run docking with single structure."""
        logger.info(f"Docking {len(dataset)} compounds to {Path(structure).name}...")
        
        results = []
        
        for idx, row in dataset.iterrows():
            compound_id = row['compound_id']
            compound = row['compound_data']
            is_active = row['is_active']
            
            try:
                # Dock compound
                docking_result = self.docking_manager.dock_compound(
                    compound=compound,
                    protein=structure,
                    grid_params=grid_params
                )
                
                score = docking_result['best_score']
                
            except Exception as e:
                logger.warning(f"Docking failed for {compound_id}: {e}")
                score = 0.0  # Assign worst score
            
            results.append({
                'compound_id': compound_id,
                'is_active': is_active,
                'docking_score': score
            })
            
            # Progress
            if (idx + 1) % 50 == 0 or (idx + 1) == len(dataset):
                logger.info(f"  Progress: {idx+1}/{len(dataset)}")
        
        results_df = pd.DataFrame(results)
        logger.info(f"Single structure docking complete")
        
        return results_df
    
    def _run_ensemble_experiment(self, dataset, structures, grid_params):
        """Run docking with ensemble."""
        logger.info(f"Docking {len(dataset)} compounds to {len(structures)} structures...")
        
        results = []
        
        for idx, row in dataset.iterrows():
            compound_id = row['compound_id']
            compound = row['compound_data']
            is_active = row['is_active']
            
            try:
                # Dock to ensemble
                ensemble_result = self.docking_manager.dock_ensemble(
                    compound=compound,
                    proteins=structures,
                    grid_params=grid_params
                )
                
                # Take best score across ensemble
                score = ensemble_result['best_score']
                best_structure = ensemble_result.get('best_structure_id', None)
                individual_scores = ensemble_result.get('individual_scores', [])
                
            except Exception as e:
                logger.warning(f"Ensemble docking failed for {compound_id}: {e}")
                score = 0.0
                best_structure = None
                individual_scores = []
            
            results.append({
                'compound_id': compound_id,
                'is_active': is_active,
                'docking_score': score,
                'best_structure': best_structure,
                'individual_scores': individual_scores
            })
            
            # Progress
            if (idx + 1) % 50 == 0 or (idx + 1) == len(dataset):
                logger.info(f"  Progress: {idx+1}/{len(dataset)}")
        
        results_df = pd.DataFrame(results)
        logger.info(f"Ensemble docking complete")
        
        return results_df
    
    def _make_recommendation(self, results):
        """Make final recommendation based on results."""
        single_auc = results['single_metrics']['AUC']
        ensemble_auc = results['ensemble_metrics']['AUC']
        
        auc_improvement = (ensemble_auc - single_auc) / single_auc * 100
        auc_significant = results['bootstrap']['auc']['significant']
        
        single_ef1 = results['single_metrics']['EF_1%']
        ensemble_ef1 = results['ensemble_metrics']['EF_1%']
        ef1_improvement = (ensemble_ef1 - single_ef1) / single_ef1 * 100 if single_ef1 > 0 else 0
        ef1_significant = results['bootstrap']['ef1']['significant']
        
        logger.info(f"\nAUC:  Single={single_auc:.3f}, Ensemble={ensemble_auc:.3f} "
                   f"({auc_improvement:+.1f}%, {'sig' if auc_significant else 'n.s.'})")
        logger.info(f"EF1%: Single={single_ef1:.1f}, Ensemble={ensemble_ef1:.1f} "
                   f"({ef1_improvement:+.1f}%, {'sig' if ef1_significant else 'n.s.'})")
        
        # Decision logic
        if ensemble_auc > single_auc and auc_significant:
            recommendation = "USE_ENSEMBLE"
            reason = "Ensemble provides statistically significant improvement in AUC."
            explanation = (
                "✅ RECOMMENDED: Use ENSEMBLE docking\n\n"
                f"   • AUC improvement: {auc_improvement:+.1f}% (p < 0.05)\n"
                f"   • EF1% improvement: {ef1_improvement:+.1f}%\n"
                "   • The 5× computational cost is JUSTIFIED\n\n"
                "   Ensemble docking captures protein flexibility better and provides\n"
                "   superior enrichment of actives over single structure approach."
            )
        elif ensemble_auc > single_auc and not auc_significant:
            recommendation = "USE_SINGLE"
            reason = "Ensemble shows improvement but not statistically significant."
            explanation = (
                "⚠️  RECOMMENDED: Use SINGLE structure (marginal benefit)\n\n"
                f"   • AUC improvement: {auc_improvement:+.1f}% (NOT significant)\n"
                f"   • EF1% improvement: {ef1_improvement:+.1f}%\n"
                "   • Save 5× computational resources\n\n"
                "   While ensemble shows slight improvement, the difference is not\n"
                "   statistically significant. Use single structure for efficiency."
            )
        else:
            recommendation = "USE_SINGLE"
            reason = "Ensemble does not improve over single structure."
            explanation = (
                "❌ RECOMMENDED: Use SINGLE structure\n\n"
                f"   • AUC difference: {auc_improvement:+.1f}%\n"
                "   • No benefit from ensemble\n"
                "   • Use faster single-structure approach\n\n"
                "   Ensemble docking does not provide improvement for this system.\n"
                "   Single structure is more efficient."
            )
        
        print("\n" + explanation)
        
        return {
            'recommendation': recommendation,
            'reason': reason,
            'auc_improvement_pct': auc_improvement,
            'auc_significant': auc_significant,
            'ef1_improvement_pct': ef1_improvement,
            'ef1_significant': ef1_significant
        }
    
    def _generate_all_plots(self, single_analysis, ensemble_analysis, bootstrap_results):
        """Generate all comparison plots."""
        logger.info("Generating plots...")
        
        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['figure.dpi'] = 150
        
        # Create figure with subplots
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        
        # 1. ROC Curves
        ax1 = fig.add_subplot(gs[0, 0])
        self._plot_roc_curves(ax1, single_analysis, ensemble_analysis)
        
        # 2. Enrichment Curves
        ax2 = fig.add_subplot(gs[0, 1])
        self._plot_enrichment_curves(ax2, single_analysis, ensemble_analysis)
        
        # 3. Metric Comparison
        ax3 = fig.add_subplot(gs[0, 2])
        self._plot_metric_comparison(ax3, single_analysis, ensemble_analysis)
        
        # 4. Bootstrap AUC Distribution
        ax4 = fig.add_subplot(gs[1, 0])
        self._plot_bootstrap_distribution(ax4, bootstrap_results, 'auc', 'AUC')
        
        # 5. Bootstrap EF1% Distribution
        ax5 = fig.add_subplot(gs[1, 1])
        self._plot_bootstrap_distribution(ax5, bootstrap_results, 'ef1', 'EF 1%')
        
        # 6. Score Distribution
        ax6 = fig.add_subplot(gs[1, 2])
        self._plot_score_distributions(ax6, single_analysis, ensemble_analysis)
        
        # 7. Early Enrichment Detail
        ax7 = fig.add_subplot(gs[2, :])
        self._plot_early_enrichment_detail(ax7, single_analysis, ensemble_analysis)
        
        # Save
        output_file = self.output_dir / 'ensemble_validation_plots.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"  Saved: {output_file}")
        
        plt.close()
    
    def _plot_roc_curves(self, ax, single_analysis, ensemble_analysis):
        """Plot ROC curves."""
        single_roc = single_analysis.get_roc_curve()
        ensemble_roc = ensemble_analysis.get_roc_curve()
        
        single_auc = single_analysis.calculate_all_metrics()['AUC']
        ensemble_auc = ensemble_analysis.calculate_all_metrics()['AUC']
        
        ax.plot(single_roc['fpr'], single_roc['tpr'], 
                label=f'Single (AUC={single_auc:.3f})', linewidth=2)
        ax.plot(ensemble_roc['fpr'], ensemble_roc['tpr'], 
                label=f'Ensemble (AUC={ensemble_auc:.3f})', linewidth=2)
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='Random')
        
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title('ROC Curves')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def _plot_enrichment_curves(self, ax, single_analysis, ensemble_analysis):
        """Plot enrichment curves."""
        single_curve = single_analysis.get_enrichment_curve()
        ensemble_curve = ensemble_analysis.get_enrichment_curve()
        
        ax.plot(single_curve['fraction_screened'] * 100,
                single_curve['enrichment_factor'],
                label='Single', linewidth=2)
        ax.plot(ensemble_curve['fraction_screened'] * 100,
                ensemble_curve['enrichment_factor'],
                label='Ensemble', linewidth=2)
        ax.axhline(y=1, color='k', linestyle='--', alpha=0.3, label='Random')
        
        ax.set_xlabel('% Database Screened')
        ax.set_ylabel('Enrichment Factor')
        ax.set_title('Enrichment Curves')
        ax.set_xlim(0, 100)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def _plot_metric_comparison(self, ax, single_analysis, ensemble_analysis):
        """Plot bar chart of metric comparison."""
        single_metrics = single_analysis.calculate_all_metrics()
        ensemble_metrics = ensemble_analysis.calculate_all_metrics()
        
        metrics = ['AUC', 'EF_1%', 'EF_5%', 'BEDROC_20']
        single_vals = [single_metrics[m] for m in metrics]
        ensemble_vals = [ensemble_metrics[m] for m in metrics]
        
        x = np.arange(len(metrics))
        width = 0.35
        
        ax.bar(x - width/2, single_vals, width, label='Single', alpha=0.8)
        ax.bar(x + width/2, ensemble_vals, width, label='Ensemble', alpha=0.8)
        
        ax.set_xticks(x)
        ax.set_xticklabels(metrics, rotation=0)
        ax.set_ylabel('Value')
        ax.set_title('Metric Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    def _plot_bootstrap_distribution(self, ax, bootstrap_results, metric_key, metric_name):
        """Plot bootstrap distribution."""
        data = bootstrap_results['raw_data']
        
        key_map = {
            'auc': ('aucs_a', 'aucs_b'),
            'ef1': ('ef1_a', 'ef1_b')
        }
        if metric_key in key_map:
            key_a, key_b = key_map[metric_key]
        else:
            key_a = f'{metric_key}s_a'
            key_b = f'{metric_key}s_b'
        vals_a = data[key_a]
        vals_b = data[key_b]
        
        ax.hist(vals_a, bins=30, alpha=0.6, label='Single', density=True)
        ax.hist(vals_b, bins=30, alpha=0.6, label='Ensemble', density=True)
        
        # Add mean lines
        ax.axvline(np.mean(vals_a), color='C0', linestyle='--', linewidth=2)
        ax.axvline(np.mean(vals_b), color='C1', linestyle='--', linewidth=2)
        
        ax.set_xlabel(metric_name)
        ax.set_ylabel('Density')
        ax.set_title(f'Bootstrap Distribution: {metric_name}')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def _plot_score_distributions(self, ax, single_analysis, ensemble_analysis):
        """Plot docking score distributions."""
        single_df = single_analysis.df
        ensemble_df = ensemble_analysis.df
        
        # Separate actives and decoys
        ax.hist(single_df[single_df['is_active']]['docking_score'],
                bins=30, alpha=0.5, label='Single - Actives', color='green')
        ax.hist(single_df[~single_df['is_active']]['docking_score'],
                bins=30, alpha=0.3, label='Single - Decoys', color='red')
        
        ax.set_xlabel('Docking Score')
        ax.set_ylabel('Count')
        ax.set_title('Score Distributions')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def _plot_early_enrichment_detail(self, ax, single_analysis, ensemble_analysis):
        """Plot detailed early enrichment (0-10%)."""
        single_curve = single_analysis.get_enrichment_curve()
        ensemble_curve = ensemble_analysis.get_enrichment_curve()
        
        # Filter to 0-10%
        mask_single = single_curve['fraction_screened'] <= 0.1
        mask_ensemble = ensemble_curve['fraction_screened'] <= 0.1
        
        ax.plot(single_curve[mask_single]['fraction_screened'] * 100,
                single_curve[mask_single]['enrichment_factor'],
                label='Single', linewidth=2, marker='o', markersize=3)
        ax.plot(ensemble_curve[mask_ensemble]['fraction_screened'] * 100,
                ensemble_curve[mask_ensemble]['enrichment_factor'],
                label='Ensemble', linewidth=2, marker='s', markersize=3)
        
        ax.axhline(y=1, color='k', linestyle='--', alpha=0.3)
        
        ax.set_xlabel('% Database Screened')
        ax.set_ylabel('Enrichment Factor')
        ax.set_title('Early Enrichment (0-10%)')
        ax.set_xlim(0, 10)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def _save_results(self, results):
        """Save results to files."""
        logger.info("Saving results...")
        
        # Save metrics as JSON
        metrics_file = self.output_dir / 'metrics.json'
        def _to_jsonable(value):
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, np.bool_):
                return bool(value)
            if isinstance(value, np.number):
                return float(value)
            if isinstance(value, dict):
                return {k: _to_jsonable(v) for k, v in value.items()}
            if isinstance(value, (list, tuple)):
                return [_to_jsonable(v) for v in value]
            return value

        metrics_dict = {
            'single': {k: _to_jsonable(v)
                      for k, v in results['single_metrics'].items()
                      if k not in ['enrichment_curve', 'roc_curve']},
            'ensemble': {k: _to_jsonable(v)
                        for k, v in results['ensemble_metrics'].items()
                        if k not in ['enrichment_curve', 'roc_curve']},
            'recommendation': _to_jsonable(results['recommendation']),
            'bootstrap': _to_jsonable(results['bootstrap'])
        }
        
        with open(metrics_file, 'w') as f:
            json.dump(metrics_dict, f, indent=2)
        logger.info(f"  Saved: {metrics_file}")
        
        # Save raw results
        self.single_results.to_csv(self.output_dir / 'single_results.csv', index=False)
        self.ensemble_results.to_csv(self.output_dir / 'ensemble_results.csv', index=False)
        logger.info(f"  Saved: single_results.csv, ensemble_results.csv")
        
        # Save summary report
        self._save_summary_report(results)
    
    def _save_summary_report(self, results):
        """Save text summary report."""
        report_file = self.output_dir / 'SUMMARY_REPORT.txt'
        
        with open(report_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("ENSEMBLE VALIDATION EXPERIMENT - SUMMARY REPORT\n")
            f.write("="*70 + "\n\n")
            
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Dataset info
            f.write("DATASET\n")
            f.write("-"*70 + "\n")
            f.write(f"Total compounds: {results['single_metrics']['n_total']}\n")
            f.write(f"Actives: {results['single_metrics']['n_actives']}\n")
            f.write(f"Decoys: {results['single_metrics']['n_decoys']}\n\n")
            
            # Metrics
            f.write("ENRICHMENT METRICS\n")
            f.write("-"*70 + "\n")
            
            single = results['single_metrics']
            ensemble = results['ensemble_metrics']
            
            f.write(f"{'Metric':<15} {'Single':>12} {'Ensemble':>12} {'Diff':>12}\n")
            f.write("-"*55 + "\n")
            
            for metric in ['AUC', 'EF_1%', 'EF_5%', 'EF_10%', 'BEDROC_20']:
                s_val = single[metric]
                e_val = ensemble[metric]
                diff = e_val - s_val
                f.write(f"{metric:<15} {s_val:>12.3f} {e_val:>12.3f} {diff:>+12.3f}\n")
            
            f.write("\n")
            
            # Statistical significance
            f.write("STATISTICAL VALIDATION (Bootstrap)\n")
            f.write("-"*70 + "\n")
            
            boot = results['bootstrap']
            f.write(f"AUC p-value: {boot['auc']['p_value']:.4f}\n")
            f.write(f"AUC significant: {boot['auc']['significant']}\n")
            f.write(f"EF1% p-value: {boot['ef1']['p_value']:.4f}\n")
            f.write(f"EF1% significant: {boot['ef1']['significant']}\n\n")
            
            # Recommendation
            f.write("RECOMMENDATION\n")
            f.write("="*70 + "\n")
            rec = results['recommendation']
            f.write(f"Decision: {rec['recommendation']}\n")
            f.write(f"Reason: {rec['reason']}\n")
            f.write(f"AUC improvement: {rec['auc_improvement_pct']:+.1f}%\n")
            f.write(f"EF1% improvement: {rec['ef1_improvement_pct']:+.1f}%\n")
        
        logger.info(f"  Saved: {report_file}")
