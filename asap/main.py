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
    args = docopt.docopt(__doc__)

    # Extract command
    command = args.pop('<command>')
    # Extract options and argument
    command_args = args.pop('<args>')

    if command == "shrake_rupley":
        import shrake_rupley
        shrake_rupley.shrake_rupley_cli(command_args)

if __name__ == '__main__':
    main()