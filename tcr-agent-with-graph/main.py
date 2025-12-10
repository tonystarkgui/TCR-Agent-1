import os
from src.agent import ImmunotherapyAgent

def main():
    config_path = "configs/solid_tumor_job.yaml"
    
    # Ensure results directory exists (handled in agent but good to have here too or just rely on agent)
    if not os.path.exists("results"):
        os.makedirs("results")
        
    agent = ImmunotherapyAgent(config_path)
    agent.run_design_cycle()

if __name__ == "__main__":
    main()
