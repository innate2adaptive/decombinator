from .decombine import decombinator
from .collapse import collapsinator
from .translate import cdr3translator
from .io import write_out_intermediate, write_out_translated, cli_args
from datetime import datetime
from typing import Optional, Any
from importlib import metadata


def run(
    args: Optional[dict[str, Any]] = None,
    cli_args: Optional[dict[str, Any]] = None,
):
    """
    Run the Decombinator pipeline
    """
    startTime = datetime.now()
    if not cli_args:
        input = args
    else:
        input = cli_args
    # Run pipline, ovewriting data after each function call to save memory
    data = decombinator(input)
    if not input["dontsave"]:
        write_out_intermediate(data, input, ".n12")
    print("Decombinator complete...")

    data = collapsinator(data=data, inputargs=input)
    if not input["dontsave"]:
        write_out_intermediate(data, input, ".freq")
    print("Collapsinator complete...")

    data = cdr3translator(data=data, inputargs=input)
    print("CDR3translator complete...")

    if not input["dontsave"]:
        write_out_translated(data, input)
    print(f"Pipeline complete in {datetime.now() - startTime}")


def main():
    input = cli_args()
    if input["command"] == "decombine":
        data = decombinator(input)
        write_out_intermediate(data, input, ".n12")
    elif input["command"] == "collapse":
        data = collapsinator(inputargs=input)
        write_out_intermediate(data, input, ".freq")
    elif input["command"] == "translate":
        data = cdr3translator(inputargs=input)
        write_out_translated(data, input)
    else:
        run(cli_args=input)


# If called from the CLI
if __name__ == "__main__":
    main()
