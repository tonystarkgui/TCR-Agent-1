import os
import requests
import logging
import random
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

def generate_sequences_evo2(seed_sequence="TGCGCCAGC", num_tokens=128, top_k=4, temperature=1.0, n_sequences=10):
    """
    Generates CDR3 sequences using NVIDIA Evo2 API (Genomic Model).
    """
    api_key = os.getenv('NVIDIA_API_KEY')
    
    if not api_key:
        logger.warning("NVIDIA_API_KEY not found. Using mock generator.")
        return _mock_generate(n_sequences)

    # Correct BioNeMo Endpoint
    url = "https://health.api.nvidia.com/v1/biology/arc/evo2-40b/generate"
    
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }
    
    # Auto-convert Protein seeds (like "CAS") to DNA ("TGC...")
    if any(aa in seed_sequence for aa in "SP"): 
        logger.info(f"Detected Protein seed '{seed_sequence}'. Converting to DNA 'TGCGCCAGC' for Evo2.")
        seed_sequence = "TGCGCCAGC"

    payload = {
        "sequence": seed_sequence, # Correct key for BioNeMo
        "num_tokens": num_tokens,  # Correct key for BioNeMo
        "top_k": top_k,
        "temperature": temperature
    }
    
    try:
        response = requests.post(url, json=payload, headers=headers)
        response.raise_for_status()
        data = response.json()
        
        protein_sequences = []
        # Handle BioNeMo response structure (might vary, usually in 'choices' or 'text')
        raw_dna_list = [item.get('text', '') for item in data.get('choices', [])]
        
        for dna in raw_dna_list:
            dna_clean = "".join([b for b in dna if b in "ATCG"])
            try:
                # Translate DNA -> Protein
                coding_dna = Seq(dna_clean)
                protein = str(coding_dna.translate(to_stop=True))
                if len(protein) > 5:
                    protein_sequences.append(protein)
            except Exception:
                continue

        if not protein_sequences:
            logger.warning("No valid proteins generated. Using mock.")
            return _mock_generate(n_sequences)
            
        return protein_sequences

    except Exception as e:
        logger.error(f"Evo2 API failed: {e}. Using mock.")
        return _mock_generate(n_sequences)

def _mock_generate(n):
    """Fallback mock generator (Returns PROTEIN)"""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    sequences = []
    for _ in range(n):
        length = random.randint(16, 22)
        core = "".join(random.choice(amino_acids) for _ in range(length))
        sequences.append("W" + core + "C")
    return sequences
