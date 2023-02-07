import argparse
import sys


def write_output(output_path: str, contents: str):
    with open(output_path, "w") as output:
        output.write(contents)


def read_input(input_path: str) -> str:
    with open(input_path, "r") as input:
        file_contents = input.read()

    return file_contents


def parse_args(sys_args: list) -> argparse.Namespace:
    parser = argparse.ArgumentParser(sys_args, description="Attempts to remove invisible Unicode chars from text content")
    parser.add_argument("input", type=str, help="Input file that we want to strip of invisible Unicode characters")
    parser.add_argument("output", type=str, help="Path to output file")

    return parser.parse_args()


def _main():
    args = parse_args(sys.argv[1:])
    input_path = args.input
    output_path = args.output

    contents = read_input(input_path)
    contents = contents.encode('ascii', errors="ignore").decode()
    write_output(output_path, contents)


if __name__ == "__main__":
    _main()