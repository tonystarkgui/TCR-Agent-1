import yaml
import logging
import pandas as pd
from src.generator import generate_sequences_evo2
from src.tme_filters import TMEFilter
from src.ranker import TCRRanker
# (Assume structure.py exists for docking)

class ImmunotherapyAgent:
    def __init__(self, config_path):
        with open(config_path) as f: self.config = yaml.safe_load(f)
        self.filter = TMEFilter(self.config)
        self.ranker = TCRRanker()
        self.logger = logging.getLogger("AI-AGENT")

    def run_design_cycle(self):
        print("--- Running Evo2 Agent ---")
        # 1. Generate
        raw = generate_sequences_evo2(n_sequences=self.config['design_parameters']['num_sequences'])
        
        # 2. Filter
        clean, _ = self.filter.apply_all(raw)
        print(f"Filtered {len(raw)} -> {len(clean)} candidates.")
        
        # 3. Rank
        ranked = self.ranker.score_sequences(clean)
        
        # 4. Save
        pd.DataFrame(ranked, columns=['CDR3', 'BERT_PLL_Score']).to_csv("results/candidates.csv", index=False)
        print("Saved to results/candidates.csv")
