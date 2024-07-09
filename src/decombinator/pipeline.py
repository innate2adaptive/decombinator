from .decombine import decombinator
from .collapse import collapsinator
from .translate import cdr3translator
from .io import write_out_intermediate, write_out_translated, cli_args
from datetime import datetime
from typing import Optional, Any


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
    write_out_intermediate(data, input, ".n12")
    print("Decombinator complete...")

    data = collapsinator(data, input)
    write_out_intermediate(data, input, ".freq")
    print("Collapsinator complete...")

    data = cdr3translator(data, input)
    print("CDR3translator complete...")

    write_out_translated(data, input)
    print(f"Pipeline complete in {datetime.now() - startTime}")


# If called from the CLI
if __name__ == "__main__":

    input = cli_args()
    run(cli_args=input)
