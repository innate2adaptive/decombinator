import argparse
import gzip
from importlib import metadata
import os
import pandas as pd
import typing


def handle_clash(
    parser: argparse.ArgumentParser,
    argument_name: str,
    help_text: str,
    shortcut: str,
    **kwargs,
):
    present = 0
    for action in parser._actions:
        if argument_name in action.option_strings:
            present += 1

    if not present:
        if argument_name == "-in" or argument_name == "--infile":
            parser.add_argument(
                "-in",
                "--infile",
                type=str,
                required=True,
                help=help_text,
            )
        else:
            parser.add_argument(
                shortcut,
                argument_name,
                help=help_text,
                **kwargs,
            )
    elif present > 1:
        parser.error(f"The {argument_name} argument can only be used once.")


def create_parser():
    parser = argparse.ArgumentParser(
        description="Decombinator: A fast and efficient tool for the analysis"
        " of T-cell receptor repertoire sequences produced by deep sequencing."
        " Include a positional argument to run a specific command."
        " Please see https://github.com/innate2adaptive/Decombinator/ for details."
    )

    # Add version information
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=metadata.version("decombinator"),
    )

    subparsers = parser.add_subparsers(
        dest="command", help="Available commands"
    )
    subparsers.required = False

    # Create parser for the "pipeline" command
    pipeline_parser = subparsers.add_parser(
        "pipeline", help="Run the entire Decombinator pipeline"
    )
    add_common_arguments(pipeline_parser)
    add_decombine_arguments(pipeline_parser)
    add_collapse_arguments(pipeline_parser)
    add_translate_arguments(pipeline_parser)

    # Create parser for the "decombine" command
    decombine_parser = subparsers.add_parser(
        "decombine", help="Decombine TCR reads"
    )
    add_common_arguments(decombine_parser)
    add_decombine_arguments(decombine_parser)

    # Create parser for the "collapse" command
    collapse_parser = subparsers.add_parser(
        "collapse", help="Collapse barcodes"
    )
    add_common_arguments(collapse_parser)
    add_collapse_arguments(collapse_parser)

    # Create parser for the "translate" command
    translate_parser = subparsers.add_parser(
        "translate", help="Translate Decombinator indexes"
    )
    add_common_arguments(translate_parser)
    add_translate_arguments(translate_parser)

    return parser


def add_common_arguments(parser: argparse.ArgumentParser):
    parser.add_argument(
        "-s",
        "--suppresssummary",
        action="store_true",
        help="Suppress the production of summary data log/file",
    )
    parser.add_argument(
        "-dz",
        "--dontgzip",
        action="store_true",
        help="Stop the output FASTQ files automatically being compressed with gzip",
    )
    parser.add_argument(
        "-dc",
        "--dontcount",
        action="store_true",
        help="Stop/Block printing the running count",
    )
    parser.add_argument(
        "-op",
        "--outpath",
        type=str,
        help="Path to output directory, writes to directory script was called in by default",
        required=False,
        default="",
    )
    parser.add_argument(
        "-c",
        "--chain",
        type=str,
        help="TCR chain (a/b/g/d)",
    )
    parser.add_argument(
        "-pf",
        "--prefix",
        type=str,
        default="dcr_",
        help='Specify the prefix of the output DCR file. Default = "dcr_"',
    )
    parser.add_argument(
        "-ds",
        "--dontsave",
        action="store_true",
        help="Don't save output files. For use when writing scripts which use the pipeline.",
    )


def add_decombine_arguments(parser: argparse.ArgumentParser):
    parser.add_argument(
        "-in",
        "--infile",
        type=str,
        required=True,
        help="Correctly demultiplexed/processed FASTQ file containing TCR reads",
    )
    parser.add_argument(
        "-br",
        "--bc_read",
        type=str,
        required=True,
        help="Which read has bar code (R1,R2). If used, ensure read selected is present in the same directory as the file specified by -in.",
    )
    parser.add_argument(
        "-dk", "--dontcheck", action="store_true", help="Skip the FASTQ check"
    )
    parser.add_argument(
        "-ex",
        "--extension",
        type=str,
        default="n12",
        help='Specify the file extension of the output DCR file. Default = "n12"',
    )
    parser.add_argument(
        "-or",
        "--orientation",
        type=str,
        default="reverse",
        help="Specify the orientation to search in (forward/reverse/both). Default = reverse",
    )
    parser.add_argument(
        "-tg",
        "--tags",
        type=str,
        default="extended",
        help="Specify which Decombinator tag set to use (extended or original). Default = extended",
    )
    parser.add_argument(
        "-sp",
        "--species",
        type=str,
        default="human",
        help="Specify which species TCR repertoire the data consists of (human or mouse). Default = human",
    )
    parser.add_argument(
        "-N",
        "--allowNs",
        action="store_true",
        help="Whether to allow VJ rearrangements containing ambiguous base calls ('N'). Default = False",
    )
    parser.add_argument(
        "-ln",
        "--lenthreshold",
        type=int,
        default=130,
        help="Acceptable threshold for inter-tag (V to J) sequence length. Default = 130",
    )
    parser.add_argument(
        "-tfdir",
        "--tagfastadir",
        type=str,
        default="Decombinator-Tags-FASTAs",
        help='Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis. Default = "Decombinator-Tags-FASTAs".',
    )
    parser.add_argument(
        "-nbc",
        "--nobarcoding",
        action="store_true",
        help="Option to run Decombinator without barcoding, i.e. so as to run on data produced by any protocol.",
    )
    parser.add_argument(
        "-bl",
        "--bclength",
        type=int,
        default=42,
        help="Length of barcode sequence, if applicable. Default is set to 42 bp.",
    )


def add_collapse_arguments(parser: argparse.ArgumentParser):
    handle_clash(
        parser,
        argument_name="--infile",
        shortcut="-in",
        help_text="File containing raw verbose Decombinator output, i.e. 5 part classifier plus barcode and inter-tag sequence and quality strings",
    )
    parser.add_argument(
        "-mq",
        "--minbcQ",
        type=int,
        default=20,
        help="Minimum quality score that barcode nucleotides should be to for that rearrangement to be retained. Default = 20.",
    )
    parser.add_argument(
        "-bm",
        "--bcQbelowmin",
        type=int,
        default=1,
        help="Number of nucleotides per barcode whose quality score are allowed to be below -mq and still be retained. Default = 1.",
    )
    parser.add_argument(
        "-aq",
        "--avgQthreshold",
        type=int,
        default=30,
        help="Average quality threshold that barcode sequences must remain above for rearrangements to be retained. Default = 30",
    )
    parser.add_argument(
        "-lv",
        "--percentlevdist",
        type=int,
        default=10,
        help="Percentage Levenshtein distance that is allowed to estimate whether two sequences within a barcode are derived from the same originator molecule. Default = 10",
    )
    parser.add_argument(
        "-bc",
        "--bcthreshold",
        type=int,
        default=2,
        help="Number of sequence edits that are allowed to consider two barcodes to be derived from same originator during clustering. Default = 2.",
    )
    handle_clash(
        parser,
        argument_name="--extension",
        shortcut="-ex",
        help_text="Specify the file extension of the output DCR file. Default = 'freq'",
        default="freq",
        type=str,
        required=False,
    )
    handle_clash(
        parser,
        argument_name="--allowNs",
        shortcut="-N",
        help_text="Used to allow VJ rearrangements containing ambiguous base calls ('N')",
        action="store_true",
    )
    handle_clash(
        parser,
        argument_name="--lenthreshold",
        shortcut="-ln",
        help_text="Acceptable threshold for inter-tag (V to J) sequence length",
        default=130,
        type=int,
        required=False,
    )
    parser.add_argument(
        "-di",
        "--dontcheckinput",
        action="store_true",
        help="Override the input file sanity check",
    )
    parser.add_argument(
        "-bd",
        "--barcodeduplication",
        action="store_true",
        help="Optionally output a file containing the final list of clustered barcodes, and their frequencies",
    )
    parser.add_argument(
        "-pb",
        "--positionalbarcodes",
        action="store_true",
        help="Instead of inferring random barcode sequences from their context relative to spacer sequences, just take the sequence at the default positions. Useful to salvage runs when R2 quality is terrible.",
    )
    parser.add_argument(
        "-ol",
        "--oligo",
        type=str,
        required=True,
        default="m13",
        help='Choose experimental oligo for correct identification of spacers ["M13", "I8", "I8_single", "NEBIO", "TAKARA"] (default: M13)',
    )
    parser.add_argument(
        "-wc",
        "--writeclusters",
        action="store_true",
        help="Write cluster data to separate cluster files",
    )
    parser.add_argument(
        "-uh",
        "--UMIhistogram",
        action="store_true",
        help="Creates histogram of average UMI cluster sizes",
    )


def add_translate_arguments(parser: argparse.ArgumentParser):
    handle_clash(
        parser,
        argument_name="--infile",
        shortcut="-in",
        help_text="File containing 5 part classifier plus barcode and inter-tag sequence and quality strings",
    )
    handle_clash(
        parser,
        argument_name="--species",
        shortcut="-sp",
        help_text="Specify which species TCR repertoire the data consists of (human or mouse). Default = human",
        default="human",
        type=str,
        required=False,
    )
    handle_clash(
        parser,
        argument_name="--tags",
        shortcut="-tg",
        help_text="Specify which Decombinator tag set to use (extended or original). Default = extended",
        default="extended",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-npf",
        "--nonproductivefilter",
        action="store_true",
        help="Filter out non-productive reads from the output",
    )
    handle_clash(
        parser,
        argument_name="--tagfastadir",
        shortcut="-tfdir",
        help_text="Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis. Default = 'Decombinator-Tags-FASTAs'",
        default="Decombinator-Tags-FASTAs",
        type=str,
        required=False,
    )
    handle_clash(
        parser,
        argument_name="--nobarcoding",
        shortcut="-nbc",
        help_text="Option to run CD3translator without barcoding, i.e. so as to run on data produced by any protocol.",
        action="store_true",
    )


def cli_args():
    parser = create_parser()
    return vars(parser.parse_args())


def create_args_dict(
    infile: str,
    chain: str,
    bc_read: str,
    suppresssummary: bool = False,
    dontgzip: bool = False,
    dontcheck: bool = False,
    dontcount: bool = False,
    extension: str = "n12",
    prefix: str = "dcr_",
    orientation: str = "reverse",
    tags: str = "extended",
    species: str = "human",
    allowNs: bool = False,
    lenthreshold: int = 130,
    tagfastadir: str = "Decombinator-Tags-FASTAs",
    nobarcoding: bool = False,
    bclength: int = 42,
    minbcQ: int = 20,
    bcQbelowmin: int = 1,
    avgQthreshold: int = 30,
    percentlevdist: int = 10,
    bcthreshold: int = 2,
    dontcheckinput: bool = False,
    barcodeduplication: bool = False,
    positionalbarcodes: bool = False,
    oligo: str = "M13",
    writeclusters: bool = False,
    UMIhistogram: bool = False,
    nonproductivefilter: bool = False,
    outpath: str = None,
    dontsave: bool = False,
    command: str = None,
) -> dict:
    """
    Creates a function argument dictionary to be used in Decombinator,
    Collapsinator, and CDR3translator.
    """

    return {
        "infile": infile,
        "chain": chain,
        "bc_read": bc_read,
        "suppresssummary": suppresssummary,
        "dontgzip": dontgzip,
        "dontcheck": dontcheck,
        "dontcount": dontcount,
        "extension": extension,
        "prefix": prefix,
        "orientation": orientation,
        "tags": tags,
        "species": species,
        "allowNs": allowNs,
        "lenthreshold": lenthreshold,
        "tagfastadir": tagfastadir,
        "nobarcoding": nobarcoding,
        "bclength": bclength,
        "minbcQ": minbcQ,
        "bcQbelowmin": bcQbelowmin,
        "avgQthreshold": avgQthreshold,
        "percentlevdist": percentlevdist,
        "bcthreshold": bcthreshold,
        "dontcheckinput": dontcheckinput,
        "barcodeduplication": barcodeduplication,
        "positionalbarcodes": positionalbarcodes,
        "oligo": oligo,
        "writeclusters": writeclusters,
        "UMIhistogram": UMIhistogram,
        "nonproductivefilter": nonproductivefilter,
        "outpath": outpath,
        "dontsave": dontsave,
        "command": command,
    }


def sort_permissions(fl):
    """
    Need to ensure proper file permissions on output data.
    If users are running pipeline through Docker might otherwise require root access
    :param fl: The file to sort permissions on
    :return: Nothing: script edits permissions where appropriate, if possible
    """

    if oct(os.stat(fl).st_mode)[4:] != "666":
        os.chmod(fl, 0o666)


def write_out_intermediate(data: list, inputargs: dict, suffix: str):
    chain = inputargs["chain"]
    chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}
    filename_id = os.path.basename(inputargs["infile"]).split(".")[0]
    if inputargs["command"] in ["collapse", "translate"]:
        outfilename = inputargs["outpath"] + f"{filename_id}" + suffix
    else:
        outfilename = (
            inputargs["outpath"]
            + inputargs["prefix"]
            + f"{filename_id}"
            + f"_{chainnams[chain.lower()]}"
            + suffix
        )
    with open(outfilename, "w") as outfile:
        for line in data:
            outfile.write(", ".join(map(str, line)) + "\n")

    if not inputargs["dontgzip"]:
        print("Compressing intermediate output file to", outfilename + ".gz")

        with (
            open(outfilename) as infile,
            gzip.open(outfilename + ".gz", "wt") as outfile,
        ):
            outfile.writelines(infile)
        os.unlink(outfilename)

        outfilenam = outfilename + ".gz"

    else:
        outfilenam = outfilename

    sort_permissions(outfilenam)


def write_out_translated(data: pd.DataFrame, inputargs: dict):
    suffix = ".tsv"
    chain = inputargs["chain"]
    chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}
    filename_id = os.path.basename(inputargs["infile"]).split(".")[0]
    if inputargs["command"] in ["collapse", "translate"]:
        outfilename = inputargs["outpath"] + f"{filename_id}" + suffix
    else:
        outfilename = (
            inputargs["outpath"]
            + inputargs["prefix"]
            + f"{filename_id}"
            + f"_{chainnams[chain.lower()]}"
            + suffix
        )
    data.to_csv(f"{outfilename}", sep="\t", index=False)

    if not inputargs["dontgzip"]:
        print("Compressing pipeline output file to", outfilename + ".gz")

        with (
            open(outfilename) as infile,
            gzip.open(outfilename + ".gz", "wt") as outfile,
        ):
            outfile.writelines(infile)
        os.unlink(outfilename)

        outfilenam = outfilename + ".gz"

    else:
        outfilenam = outfilename

    sort_permissions(outfilenam)
