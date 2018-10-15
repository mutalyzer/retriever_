"""
CLI entry point.
"""

import argparse

from . import usage, version


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        description=usage[0], epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', action='version', version=version(parser.prog))

    args = parser.parse_args()