from Bio.SeqUtils.ProtParam import ProteinAnalysis

class TMEFilter:
    def __init__(self, config):
        self.constraints = config.get('tme_constraints', [])

    def check_exhaustion_risk(self, sequence):
        """
        High hydrophobicity in CDR3 often leads to 'Tonic Signaling',
        driving T-cells to exhaustion even without antigen.
        We filter out 'sticky' sequences.
        """
        analysis = ProteinAnalysis(sequence)
        gravy_score = analysis.gravy() # Grand Average of Hydropathy
        
        # Threshold: Lower is less hydrophobic (better)
        # If score > threshold, risk of tonic signaling is HIGH.
        threshold = next((c['threshold'] for c in self.constraints if c['constraint'] == 'tonic_signaling_risk'), 0.5)
        
        if gravy_score > threshold:
            return False, f"Rejected: High Tonic Signaling Risk (GRAVY: {gravy_score:.2f})"
        return True, "Passed"

    def check_stability(self, sequence):
        """
        Checks for instability motifs (e.g., NG deamidation sites) 
        that would degrade in the acidic TME.
        """
        banned = next((c['motifs'] for c in self.constraints if c['constraint'] == 'motif_ban'), [])
        for motif in banned:
            if motif in sequence:
                return False, f"Rejected: Unstable Motif ({motif})"
        return True, "Passed"

    def apply_all(self, sequences):
        passed_sequences = []
        logs = []
        
        for seq in sequences:
            # 1. Exhaustion Check
            ok_ex, msg_ex = self.check_exhaustion_risk(seq)
            if not ok_ex:
                logs.append(msg_ex)
                continue
                
            # 2. Stability Check
            ok_stab, msg_stab = self.check_stability(seq)
            if not ok_stab:
                logs.append(msg_stab)
                continue
                
            passed_sequences.append(seq)
            
        return passed_sequences, logs
