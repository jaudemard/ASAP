"""
usage:
    asap <command> [<args>...]

options:
    -h --help  Show help
commands:
    shrake_rupley   Determine Solvent Accessible Surface Area with a Shrake and Rupley based method.
"""
"""
usage:
    asap <command> [<args>...]

options:
    -h --help  Show help
commands:
    shrake_rupley   Determine Solvent Accessible Surface Area with a Shrake and Rupley based method.
"""
import docopt
import importlib
import sys

def main(args=None):
    # Parse the command-line arguments
    args = docopt.docopt(__doc__)

    # Extract command and arguments
    command = args.pop('<command>')
    command_args = args.pop('<args>')
    print(command_args)

    try:
        # Import the command module
        command_module = importlib.import_module(command)
    except ModuleNotFoundError:
        print(f"Unrecognized command '{command}'")
        print(__doc__)  # Print the usage guide
        sys.exit(1)

    # Build the CLI function name
    cli_function_name = f"{command}_cli"
    
    # Check if the CLI function exist
    if hasattr(command_module, cli_function_name):
        cli_function = getattr(command_module, cli_function_name)
        if ("-h" or "--help") in command_args:
            command_module.command_help()
        else:
            # Execute the CLI function with user's arguments
            cli_function(command_args)
    else:
        print(f"'{command}' does not have a CLI interface.")
        print(__doc__)
        sys.exit(1)

if __name__ == '__main__':
    main()