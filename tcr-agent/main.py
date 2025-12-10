import argparse
import sys
import os
from src.agent import ImmunotherapyAgent

def main():
    parser = argparse.ArgumentParser(description="TCR Agent Runner")
    parser.add_argument("--config", type=str, required=True, help="Path to config yaml")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: Config file not found at {args.config}")
        sys.exit(1)

    try:
        agent = ImmunotherapyAgent(args.config)
        agent.run_design_cycle()
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Execution failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
