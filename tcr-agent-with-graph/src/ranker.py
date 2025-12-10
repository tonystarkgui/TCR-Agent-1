import torch
from transformers import BertForMaskedLM, BertTokenizer
import logging
import numpy as np

logger = logging.getLogger(__name__)

class TCRRanker:
    def __init__(self, model_path="wukevin/tcr-bert"):
        self.model_path = model_path
        self.tokenizer = None
        self.model = None
        self._load_model()

    def _load_model(self):
        try:
            logger.info(f"Loading TCR-BERT model from {self.model_path}...")
            self.tokenizer = BertTokenizer.from_pretrained(self.model_path)
            self.model = BertForMaskedLM.from_pretrained(self.model_path)
            self.model.eval()
        except Exception as e:
            logger.error(f"Failed to load TCR-BERT model: {e}. Ranking will be mocked.")
            self.model = None

    def score_sequences(self, sequences):
        """
        Scores sequences using Pseudo-Log-Likelihood (PLL).
        Higher PLL (closer to 0) is better.
        """
        if not self.model:
            return self._mock_score(sequences)

        scored_seqs = []
        for seq in sequences:
            pll = self._calculate_pll(seq)
            scored_seqs.append((seq, pll))
        
        # Sort by PLL descending (higher is better)
        scored_seqs.sort(key=lambda x: x[1], reverse=True)
        return scored_seqs

    def _calculate_pll(self, sequence):
        """
        Calculates Pseudo-Log-Likelihood for a single sequence.
        """
        # Tokenize
        # TCR-BERT expects spaced amino acids: "C A S ..."
        spaced_seq = " ".join(list(sequence))
        inputs = self.tokenizer(spaced_seq, return_tensors="pt")
        input_ids = inputs["input_ids"]
        labels = input_ids.clone()

        total_log_prob = 0.0
        length = len(sequence) # Effective length (ignoring CLS/SEP for averaging if desired, but here we just sum/avg)
        
        # Loop through each token (skipping CLS and SEP)
        # input_ids[0] is [CLS, A, A, ..., SEP]
        # range(1, len - 1)
        
        seq_len = input_ids.size(1)
        
        with torch.no_grad():
            for i in range(1, seq_len - 1): # Skip CLS (0) and SEP (last)
                # Mask the current token
                orig_token_id = input_ids[0, i].item()
                input_ids[0, i] = self.tokenizer.mask_token_id
                
                # Predict
                outputs = self.model(input_ids)
                predictions = outputs.logits
                
                # Get log prob of original token
                # predictions shape: [1, seq_len, vocab_size]
                token_logits = predictions[0, i, :]
                token_log_probs = torch.nn.functional.log_softmax(token_logits, dim=0)
                token_pll = token_log_probs[orig_token_id].item()
                
                total_log_prob += token_pll
                
                # Restore token
                input_ids[0, i] = orig_token_id
        
        # Average PLL
        avg_pll = total_log_prob / length if length > 0 else -999.0
        return avg_pll

    def _mock_score(self, sequences):
        logger.warning("Using mock ranking.")
        # Return random scores
        import random
        scored = [(seq, -1.0 * random.random()) for seq in sequences]
        scored.sort(key=lambda x: x[1], reverse=True)
        return scored
