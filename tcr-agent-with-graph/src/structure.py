import os
import logging

logger = logging.getLogger(__name__)

def prepare_docking_job(sequences, target_name):
    """
    Prepares input for TCRDock for the top sequences.
    """
    logger.info(f"Preparing TCRDock job for {len(sequences)} sequences against {target_name}...")
    
    # Create a job directory
    job_dir = "results/docking_jobs"
    os.makedirs(job_dir, exist_ok=True)
    
    job_file = os.path.join(job_dir, f"{target_name}_job.csv")
    
    with open(job_file, "w") as f:
        f.write("id,cdr3_sequence,target\n")
        for i, (seq, score) in enumerate(sequences):
            f.write(f"candidate_{i+1},{seq},{target_name}\n")
            
    logger.info(f"Docking job file created at: {job_file}")
    
    # Generate a mock command string to run TCRDock
    cmd = f"tcrdock --input {job_file} --output results/docking_output"
    logger.info(f"Run command: {cmd}")
    
    return job_file
