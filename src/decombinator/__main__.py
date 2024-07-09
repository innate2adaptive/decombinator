from .io import cli_args
from .pipeline import run

input = cli_args()
run(cli_args=input)
