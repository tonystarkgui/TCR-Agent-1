import yaml
import logging
import os
import pandas as pd
from src.generator import generate_sequences_evo2
from src.tme_filters import TMEFilter
from src.ranker import TCRRanker
from src.structure import prepare_docking_job

class ImmunotherapyAgent:
    def __init__(self, config_path):
        with open(config_path) as f:
            self.config = yaml.safe_load(f)
        self.logger = self._setup_logger()
        self.filter_module = TMEFilter(self.config)
        
        # Initialize Ranker
        bert_path = self.config.get('bert_model_path', 'wukevin/tcr-bert')
        self.ranker = TCRRanker(bert_path)

    def _setup_logger(self):
        logging.basicConfig(level=logging.INFO, format='[AI-AGENT] %(message)s')
        return logging.getLogger(__name__)

    def run_design_cycle(self):
        target = self.config['target']['name']
        self.logger.info(f"--- INIT: Designing T-cells for Solid Tumor Target: {target} ---")
        
        # 1. Generation (NVIDIA Evo2)
        self.logger.info("Step 1: Generating candidates via NVIDIA Evo2...")
        evo_params = self.config.get('evo2_parameters', {})
        n_seqs = self.config['design_parameters']['num_sequences']
        
        raw_seqs = generate_sequences_evo2(
            num_tokens=20, # Default
            top_k=evo_params.get('top_k', 4),
            temperature=evo_params.get('temperature', 1.0),
            n_sequences=n_seqs
        )
        self.logger.info(f"Generated {len(raw_seqs)} raw CDR3 sequences.")

        # 2. TME/Exhaustion Filtering (CRITICAL)
        self.logger.info("Step 2: Filtering for TME Survival & Low Exhaustion Risk...")
        clean_seqs, logs = self.filter_module.apply_all(raw_seqs)
        
        rejection_rate = (1 - len(clean_seqs)/len(raw_seqs)) * 100 if raw_seqs else 0
        self.logger.info(f"TME Filter removed {rejection_rate:.1f}% of candidates (Risk of Tonic Signaling/Instability).")
        self.logger.info(f"Candidates remaining: {len(clean_seqs)}")

        if not clean_seqs:
            self.logger.warning("No sequences survived TME filtering. Aborting.")
            return

        # 3. Ranking (TCR-BERT)
        self.logger.info("Step 3: Ranking candidates with TCR-BERT (PLL Scoring)...")
        ranked_candidates = self.ranker.score_sequences(clean_seqs)
        
        # Log top 3
        self.logger.info("Top 3 Candidates by BERT Score:")
        for i, (seq, score) in enumerate(ranked_candidates[:3]):
            self.logger.info(f"  #{i+1}: {seq} (PLL: {score:.4f})")

        # 4. Structure (TCRDock Prep)
        self.logger.info("Step 4: Preparing Top 5 Candidates for TCRDock...")
        top_5 = ranked_candidates[:5]
        prepare_docking_job(top_5, target)

        # Export all results
        self.export_results(ranked_candidates)

    def export_results(self, ranked_seqs):
        # Save to CSV
        os.makedirs("results", exist_ok=True)
        
        df = pd.DataFrame(ranked_seqs, columns=['CDR3', 'BERT_PLL_Score'])
        df['Status'] = 'Ranked'
        df.to_csv("results/candidates.csv", index=False)
        self.logger.info("Saved ranked candidates to results/candidates.csv")
