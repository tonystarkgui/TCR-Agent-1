from Bio.SeqUtils.ProtParam import ProteinAnalysis

class TMEFilter:
    def __init__(self, config):
        self.constraints = config.get('tme_constraints', [])

    def apply_all(self, sequences):
        passed_sequences = []
        logs = []
        
        # Threshold defaults to 0.5 if not in config
        threshold = next((c['threshold'] for c in self.constraints if c['constraint'] == 'tonic_signaling_risk'), 0.5)

        for seq in sequences:
            # 1. Hydrophobicity Check (Tonic Signaling)
            try:
                analysis = ProteinAnalysis(seq)
                gravy = analysis.gravy()
                if gravy > threshold:
                    logs.append(f"Rejected {seq}: High GRAVY ({gravy:.2f})")
                    continue
            except:
                continue # Skip invalid sequences

            # 2. Stability Check (Motifs)
            if "NG" in seq or "NS" in seq:
                 logs.append(f"Rejected {seq}: Unstable Motif")
                 continue
                
            passed_sequences.append(seq)
            
        return passed_sequences, logs
