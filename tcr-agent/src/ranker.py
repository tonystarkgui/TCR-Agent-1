import torch
from transformers import BertForMaskedLM, BertTokenizer
import logging

logger = logging.getLogger(__name__)

class TCRRanker:
    def __init__(self, model_path="wukevin/tcr-bert-mlm-only"):
        self.tokenizer = BertTokenizer.from_pretrained(model_path)
        self.model = BertForMaskedLM.from_pretrained(model_path)
        self.model.eval()

    def score_sequences(self, sequences):
        scored = []
        for seq in sequences:
            try:
                pll = self._calculate_pll(seq)
                scored.append((seq, pll))
            except:
                scored.append((seq, -99.9))
        
        # Sort descending (Higher/Closer to 0 is better)
        scored.sort(key=lambda x: x[1], reverse=True)
        return scored

    def _calculate_pll(self, sequence):
        # Format for TCR-BERT: "C A S S ..."
        text = " ".join(list(sequence))
        inputs = self.tokenizer(text, return_tensors="pt")
        input_ids = inputs["input_ids"]
        
        total_log_prob = 0.0
        length = len(sequence)
        
        with torch.no_grad():
            for i in range(1, input_ids.size(1) - 1): # Skip CLS/SEP
                orig_id = input_ids[0, i].item()
                input_ids[0, i] = self.tokenizer.mask_token_id
                
                outputs = self.model(input_ids)
                log_probs = torch.nn.functional.log_softmax(outputs.logits[0, i], dim=0)
                total_log_prob += log_probs[orig_id].item()
                
                input_ids[0, i] = orig_id # Restore
                
        return total_log_prob / length
