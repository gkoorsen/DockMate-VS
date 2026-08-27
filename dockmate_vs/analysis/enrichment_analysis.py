"""
Enrichment analysis module for docking validation.
Calculates AUC, EF, BEDROC, and other metrics.
"""

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
from scipy import stats
from loguru import logger
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


class EnrichmentAnalysis:
    """
    Calculate enrichment metrics for docking validation.
    
    Based on best practices from:
    - Truchon & Bayly (2007) - BEDROC
    - Huang et al. (2006) - Benchmarking
    - Braga & Andrade (2013) - Metrics
    """
    
    def __init__(self, results_df):
        """
        Initialize with docking results.
        
        Parameters
        ----------
        results_df : pd.DataFrame
            Must contain columns:
            - compound_id: unique identifier
            - is_active: boolean (True for actives, False for decoys)
            - docking_score: float (lower is better for binding)
        """
        required_cols = ['compound_id', 'is_active', 'docking_score']
        missing = [col for col in required_cols if col not in results_df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        
        self.df = results_df.copy()
        
        # Sort by docking score (lower/more negative is better)
        self.df = self.df.sort_values('docking_score', ascending=True).reset_index(drop=True)
        
        self.n_total = len(self.df)
        self.n_actives = self.df['is_active'].sum()
        self.n_decoys = self.n_total - self.n_actives
        
        logger.info(f"Enrichment analysis initialized:")
        logger.info(f"  Total compounds: {self.n_total}")
        logger.info(f"  Actives: {self.n_actives}")
        logger.info(f"  Decoys: {self.n_decoys}")
        logger.info(f"  Ratio: 1:{self.n_decoys/self.n_actives:.1f}")
    
    def calculate_all_metrics(self):
        """Calculate all enrichment metrics."""
        if self.n_actives == 0:
            logger.warning("No actives in dataset - cannot calculate enrichment")
            return None
        
        metrics = {
            'AUC': self.calculate_auc(),
            'EF_0.5%': self.calculate_ef(fraction=0.005),
            'EF_1%': self.calculate_ef(fraction=0.01),
            'EF_2%': self.calculate_ef(fraction=0.02),
            'EF_5%': self.calculate_ef(fraction=0.05),
            'EF_10%': self.calculate_ef(fraction=0.10),
            'BEDROC_20': self.calculate_bedroc(alpha=20.0),
            'BEDROC_160.9': self.calculate_bedroc(alpha=160.9),
            'RIE_20': self.calculate_rie(alpha=20.0),
            'pAUC_5%': self.calculate_pauc(fraction=0.05),
            'n_actives': self.n_actives,
            'n_decoys': self.n_decoys,
            'n_total': self.n_total
        }
        
        # Add enrichment curve data
        metrics['enrichment_curve'] = self.get_enrichment_curve()
        metrics['roc_curve'] = self.get_roc_curve()
        
        return metrics
    
    def calculate_auc(self):
        """
        Calculate ROC AUC (Area Under Curve).
        
        Returns
        -------
        float : AUC value between 0 and 1 (0.5 = random, 1.0 = perfect)
        """
        if self.n_actives == 0 or self.n_decoys == 0:
            return None
        
        y_true = self.df['is_active'].astype(int).values
        # Negate scores because sklearn expects higher values for positives
        y_score = -self.df['docking_score'].values
        
        auc = roc_auc_score(y_true, y_score)
        return auc
    
    def calculate_ef(self, fraction=0.01):
        """
        Calculate Enrichment Factor at given fraction.
        
        EF_x% = (actives_in_top_x% / total_actives) / fraction
        
        Parameters
        ----------
        fraction : float
            Fraction of dataset to consider (e.g., 0.01 for 1%)
        
        Returns
        -------
        float : Enrichment factor (max value = 1/fraction)
        """
        if self.n_actives == 0:
            return 0.0
        
        # Number of compounds in top fraction
        n_top = max(1, int(self.n_total * fraction))
        
        # Get top compounds
        top_compounds = self.df.head(n_top)
        n_actives_in_top = top_compounds['is_active'].sum()
        
        # Calculate EF
        ef = (n_actives_in_top / self.n_actives) / fraction
        
        return ef
    
    def calculate_bedroc(self, alpha=20.0):
        """
        Calculate BEDROC (Boltzmann-Enhanced Discrimination of ROC).
        
        Emphasizes early recognition more than AUC.
        Reference: Truchon & Bayly (2007) J. Chem. Inf. Model.
        
        Parameters
        ----------
        alpha : float
            Controls degree of early recognition emphasis.
            alpha = 20.0 corresponds to ~8% of dataset
            alpha = 160.9 corresponds to ~1% of dataset
        
        Returns
        -------
        float : BEDROC value between 0 and 1
        """
        if self.n_actives == 0:
            return 0.0
        
        n = self.n_total
        n_actives = self.n_actives
        
        # Get ranks of actives (1-indexed)
        active_indices = np.where(self.df['is_active'])[0]
        active_ranks = active_indices + 1
        
        # Calculate RIE (Robust Initial Enhancement)
        sum_exp = np.sum(np.exp(-alpha * active_ranks / n))
        
        # Calculate normalization terms
        ra = n_actives / n
        rand_sum = ra * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
        
        rie = sum_exp / n_actives
        
        # Calculate BEDROC
        bedroc = (rie - rand_sum) / (1 - rand_sum)
        
        # Normalize to [0, 1]
        bedroc_normalized = bedroc * (1 / (1 - np.exp(-alpha)))
        
        return max(0.0, min(1.0, bedroc_normalized))
    
    def calculate_rie(self, alpha=20.0):
        """
        Calculate RIE (Robust Initial Enhancement).
        
        Similar to BEDROC but without final normalization.
        
        Parameters
        ----------
        alpha : float
            Controls degree of early recognition emphasis
        
        Returns
        -------
        float : RIE value
        """
        if self.n_actives == 0:
            return 0.0
        
        n = self.n_total
        n_actives = self.n_actives
        
        # Get ranks of actives (1-indexed)
        active_indices = np.where(self.df['is_active'])[0]
        active_ranks = active_indices + 1
        
        # Calculate RIE
        sum_exp = np.sum(np.exp(-alpha * active_ranks / n))
        rie = sum_exp / n_actives
        
        # Calculate normalization
        ra = n_actives / n
        rand_sum = ra * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
        
        rie_normalized = (rie - rand_sum) / (1 - rand_sum)
        
        return rie_normalized
    
    def calculate_pauc(self, fraction=0.05):
        """
        Calculate partial AUC (pAUC) up to given fraction.
        
        Focuses on early enrichment region of ROC curve.
        
        Parameters
        ----------
        fraction : float
            Fraction of FPR to consider (e.g., 0.05 for 5%)
        
        Returns
        -------
        float : Normalized pAUC (0.5 = random, 1.0 = perfect)
        """
        if self.n_actives == 0 or self.n_decoys == 0:
            return None
        
        y_true = self.df['is_active'].astype(int).values
        y_score = -self.df['docking_score'].values
        
        # Calculate full ROC curve
        fpr, tpr, _ = roc_curve(y_true, y_score)
        
        # Find indices up to fraction
        idx = np.where(fpr <= fraction)[0]
        
        if len(idx) == 0:
            return 0.5
        
        # Calculate partial AUC using trapezoidal rule
        fpr_partial = fpr[idx]
        tpr_partial = tpr[idx]
        
        # Add endpoint if needed
        if fpr_partial[-1] < fraction:
            # Interpolate to get TPR at exact fraction
            tpr_at_fraction = np.interp(fraction, fpr, tpr)
            fpr_partial = np.append(fpr_partial, fraction)
            tpr_partial = np.append(tpr_partial, tpr_at_fraction)
        
        # NumPy 2.x provides trapezoid; retain trapz support for older releases.
        _trap = getattr(np, "trapezoid", None) or np.trapz
        pauc = _trap(tpr_partial, fpr_partial)
        
        # Normalize: random = 0.5, perfect = 1.0
        pauc_max = fraction  # Max possible pAUC
        pauc_min = fraction * fraction / 2  # Random
        pauc_normalized = 0.5 + (pauc - pauc_min) / (pauc_max - pauc_min) * 0.5
        
        return pauc_normalized
    
    def get_enrichment_curve(self):
        """
        Get enrichment curve data for plotting.
        
        Returns
        -------
        pd.DataFrame with columns:
            - fraction_screened: fraction of database
            - enrichment_factor: EF at that fraction
            - cumulative_actives: number of actives found
        """
        cumulative_actives = np.cumsum(self.df['is_active'].values)
        fractions = np.arange(1, self.n_total + 1) / self.n_total
        
        # Enrichment factor at each point
        ef = (cumulative_actives / self.n_actives) / fractions
        
        curve_df = pd.DataFrame({
            'fraction_screened': fractions,
            'enrichment_factor': ef,
            'cumulative_actives': cumulative_actives,
            'cumulative_fraction': cumulative_actives / self.n_actives
        })
        
        return curve_df
    
    def get_roc_curve(self):
        """
        Get ROC curve data for plotting.
        
        Returns
        -------
        pd.DataFrame with columns:
            - fpr: False Positive Rate
            - tpr: True Positive Rate
            - threshold: Docking score threshold
        """
        if self.n_actives == 0 or self.n_decoys == 0:
            return None
        
        y_true = self.df['is_active'].astype(int).values
        y_score = -self.df['docking_score'].values
        
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        
        return pd.DataFrame({
            'fpr': fpr,
            'tpr': tpr,
            'threshold': thresholds
        })
    
    def compare_with(self, other_analysis, method_names=('Method A', 'Method B')):
        """
        Compare this analysis with another.
        
        Parameters
        ----------
        other_analysis : EnrichmentAnalysis
            Another analysis to compare with
        method_names : tuple
            Names for the two methods
        
        Returns
        -------
        pd.DataFrame : Comparison table
        """
        metrics_self = self.calculate_all_metrics()
        metrics_other = other_analysis.calculate_all_metrics()
        
        comparison = {
            'Metric': [],
            method_names[0]: [],
            method_names[1]: [],
            'Difference': [],
            'Percent_Change': []
        }
        
        metric_keys = ['AUC', 'EF_1%', 'EF_5%', 'EF_10%', 'BEDROC_20']
        
        for metric in metric_keys:
            val_self = metrics_self.get(metric)
            val_other = metrics_other.get(metric)
            
            if val_self is None or val_other is None:
                continue
            
            diff = val_self - val_other
            pct_change = (diff / val_other * 100) if val_other != 0 else 0
            
            comparison['Metric'].append(metric)
            comparison[method_names[0]].append(f"{val_self:.3f}")
            comparison[method_names[1]].append(f"{val_other:.3f}")
            comparison['Difference'].append(f"{diff:+.3f}")
            comparison['Percent_Change'].append(f"{pct_change:+.1f}%")
        
        return pd.DataFrame(comparison)
    
    def print_summary(self, method_name="Method"):
        """Print summary of enrichment metrics."""
        metrics = self.calculate_all_metrics()
        
        print(f"\n{'='*60}")
        print(f"ENRICHMENT SUMMARY - {method_name}")
        print(f"{'='*60}")
        print(f"Dataset: {self.n_actives} actives, {self.n_decoys} decoys")
        print(f"\nPrimary Metrics:")
        print(f"  AUC:              {metrics['AUC']:.3f}")
        print(f"  EF 1%:            {metrics['EF_1%']:.1f}")
        print(f"  EF 5%:            {metrics['EF_5%']:.1f}")
        print(f"  BEDROC (α=20):    {metrics['BEDROC_20']:.3f}")
        print(f"\nAdditional Metrics:")
        print(f"  EF 0.5%:          {metrics['EF_0.5%']:.1f}")
        print(f"  EF 2%:            {metrics['EF_2%']:.1f}")
        print(f"  EF 10%:           {metrics['EF_10%']:.1f}")
        print(f"  BEDROC (α=160.9): {metrics['BEDROC_160.9']:.3f}")
        print(f"  pAUC 5%:          {metrics['pAUC_5%']:.3f}")
        print(f"{'='*60}\n")


class BootstrapValidation:
    """Statistical validation using bootstrap resampling."""
    
    @staticmethod
    def compare_methods(results_a, results_b, n_iterations=1000, seed=42):
        """
        Compare two docking methods using bootstrap resampling.
        
        Parameters
        ----------
        results_a, results_b : pd.DataFrame
            Results from two methods (same format as EnrichmentAnalysis)
        n_iterations : int
            Number of bootstrap iterations
        seed : int
            Random seed for reproducibility
        
        Returns
        -------
        dict : Bootstrap statistics
        """
        np.random.seed(seed)
        
        n_samples = len(results_a)
        if len(results_b) != n_samples:
            raise ValueError("Both result sets must have same number of compounds")
        
        logger.info(f"Running bootstrap validation ({n_iterations} iterations)...")
        
        # Storage for bootstrap samples
        aucs_a = []
        aucs_b = []
        ef1_a = []
        ef1_b = []
        
        for i in range(n_iterations):
            # Resample with replacement
            indices = np.random.choice(n_samples, size=n_samples, replace=True)
            
            sample_a = results_a.iloc[indices].reset_index(drop=True)
            sample_b = results_b.iloc[indices].reset_index(drop=True)
            
            # Calculate metrics
            try:
                analysis_a = EnrichmentAnalysis(sample_a)
                analysis_b = EnrichmentAnalysis(sample_b)
                
                metrics_a = analysis_a.calculate_all_metrics()
                metrics_b = analysis_b.calculate_all_metrics()
                
                aucs_a.append(metrics_a['AUC'])
                aucs_b.append(metrics_b['AUC'])
                ef1_a.append(metrics_a['EF_1%'])
                ef1_b.append(metrics_b['EF_1%'])
            except Exception as e:
                logger.warning(f"Bootstrap iteration {i} failed: {e}")
                continue
            
            if (i + 1) % 100 == 0:
                logger.info(f"  Progress: {i+1}/{n_iterations}")
        
        # Calculate statistics
        aucs_a = np.array(aucs_a)
        aucs_b = np.array(aucs_b)
        ef1_a = np.array(ef1_a)
        ef1_b = np.array(ef1_b)
        
        # Confidence intervals (95%)
        auc_a_ci = np.percentile(aucs_a, [2.5, 97.5])
        auc_b_ci = np.percentile(aucs_b, [2.5, 97.5])
        ef1_a_ci = np.percentile(ef1_a, [2.5, 97.5])
        ef1_b_ci = np.percentile(ef1_b, [2.5, 97.5])
        
        # Differences
        auc_diff = aucs_b - aucs_a
        ef1_diff = ef1_b - ef1_a
        
        # P-values (one-sided: B > A)
        p_value_auc = np.sum(auc_diff <= 0) / len(auc_diff)
        p_value_ef1 = np.sum(ef1_diff <= 0) / len(ef1_diff)
        
        results = {
            'auc': {
                'method_a_mean': np.mean(aucs_a),
                'method_a_ci': auc_a_ci,
                'method_b_mean': np.mean(aucs_b),
                'method_b_ci': auc_b_ci,
                'difference_mean': np.mean(auc_diff),
                'difference_ci': np.percentile(auc_diff, [2.5, 97.5]),
                'p_value': p_value_auc,
                'significant': p_value_auc < 0.05
            },
            'ef1': {
                'method_a_mean': np.mean(ef1_a),
                'method_a_ci': ef1_a_ci,
                'method_b_mean': np.mean(ef1_b),
                'method_b_ci': ef1_b_ci,
                'difference_mean': np.mean(ef1_diff),
                'difference_ci': np.percentile(ef1_diff, [2.5, 97.5]),
                'p_value': p_value_ef1,
                'significant': p_value_ef1 < 0.05
            },
            'raw_data': {
                'aucs_a': aucs_a,
                'aucs_b': aucs_b,
                'ef1_a': ef1_a,
                'ef1_b': ef1_b
            }
        }
        
        return results
    
    @staticmethod
    def print_bootstrap_results(results, method_names=('Method A', 'Method B')):
        """Print formatted bootstrap results."""
        print(f"\n{'='*60}")
        print("BOOTSTRAP VALIDATION RESULTS")
        print(f"{'='*60}")
        
        # AUC
        auc_data = results['auc']
        print(f"\nAUC:")
        print(f"  {method_names[0]}: {auc_data['method_a_mean']:.3f} "
              f"(95% CI: {auc_data['method_a_ci'][0]:.3f}-{auc_data['method_a_ci'][1]:.3f})")
        print(f"  {method_names[1]}: {auc_data['method_b_mean']:.3f} "
              f"(95% CI: {auc_data['method_b_ci'][0]:.3f}-{auc_data['method_b_ci'][1]:.3f})")
        print(f"  Difference: {auc_data['difference_mean']:+.3f} "
              f"(95% CI: {auc_data['difference_ci'][0]:+.3f} to {auc_data['difference_ci'][1]:+.3f})")
        print(f"  P-value: {auc_data['p_value']:.4f}")
        print(f"  Significant: {'YES ✓' if auc_data['significant'] else 'NO ✗'}")
        
        # EF1%
        ef1_data = results['ef1']
        print(f"\nEF 1%:")
        print(f"  {method_names[0]}: {ef1_data['method_a_mean']:.1f} "
              f"(95% CI: {ef1_data['method_a_ci'][0]:.1f}-{ef1_data['method_a_ci'][1]:.1f})")
        print(f"  {method_names[1]}: {ef1_data['method_b_mean']:.1f} "
              f"(95% CI: {ef1_data['method_b_ci'][0]:.1f}-{ef1_data['method_b_ci'][1]:.1f})")
        print(f"  Difference: {ef1_data['difference_mean']:+.1f} "
              f"(95% CI: {ef1_data['difference_ci'][0]:+.1f} to {ef1_data['difference_ci'][1]:+.1f})")
        print(f"  P-value: {ef1_data['p_value']:.4f}")
        print(f"  Significant: {'YES ✓' if ef1_data['significant'] else 'NO ✗'}")
        
        print(f"{'='*60}\n")
