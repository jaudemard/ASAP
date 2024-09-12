import argparse
import importlib
import asap.shrake_rupley
import sys

def main():
    parser = argparse.ArgumentParser(description='Compute Solvent Accessible Surface Area of a Protein.')

    # Define the command-line arguments
    parser.add_argument('command', type=str, help='Algorithm to compute SASA: shrake_rupley')
    parser.add_argument('--pdb', type=str, help='Path to the PDB structure')
    parser.add_argument('--probe', type=float, help='Size of the probe')
    parser.add_argument('--output', type=str, default='.', help='Output directory for the output files')
    parser.add_argument('--model', type=int, default=0, help='Model used if the pdb includes several of them')
    parser.add_argument('--point', type=int, default=100, help='Number of points in the sphere lattice')
    
    args = parser.parse_args()

    # Extract command and arguments
    command_args = [
        '--pdb', str(args.pdb),
        '--probe', str(args.probe),
        '--output', str(args.output),
        '--model', str(args.model),
        '--point', str(args.point)
    ]
    # try:
    #     command_module = importlib.import_module(f"asap.shrake_rupley")
    # except ModuleNotFoundError:
    #     print(f"Unrecognized command '{command}'")
    #     parser.print_help()
    #     sys.exit(1)

    # # Build the CLI function name
    # cli_function_name = f"{command}_cli"
    # if not hasattr(command_module, cli_function_name):
    #     print(f"No CLI function '{cli_function_name}' found in module '{command}'")
    #     sys.exit(1)

    # cli_function = getattr(command_module, cli_function_name)

    # # Execute the CLI function with command_args
    asap.shrake_rupley.shrake_rupley_cli(command_args)


if __name__ == '__main__':
    main()
