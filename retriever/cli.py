"""
CLI entry point.
"""

import argparse

from . import usage, version
from .retriever import retrieve
from .ncbi import link_reference
from .parser import parse


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        description=usage[0], epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', action='version', version=version(parser.prog))

    parser.add_argument('reference', help='the reference id')

    parser.add_argument("--sizeoff", help="do not consider file size",
                        action="store_true")

    parser.add_argument("--link", help="link protein to transcript",
                        action="store_true")

    parser.add_argument("--parse", help="parse reference content",
                        action="store_true")

    args = parser.parse_args()

    if args.sizeoff:
        content, reference_type = retrieve(args.reference, not args.sizeoff)
    else:
        content, reference_type = retrieve(args.reference)

    if args.link:
        print(link_reference(args.reference))
        return

    if args.parse:
        if content:
            reference_model = parse(content, reference_type)
            print(reference_model)
        else:
            print('No content, parsing not performed.')
    else:
        print(content)
