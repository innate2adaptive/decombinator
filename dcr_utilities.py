import os
import gzip
import argparse
import pandas as pd

def cli_args():
    """args(): Obtains command line arguments which dictate the script's behaviour"""

    # Help flag
    parser = argparse.ArgumentParser(
        description='Decombinator v4.2.0: find rearranged TCR sequences in HTS data. Please go to https://innate2adaptive.github.io/Decombinator/ for more details.')
    # Decombinator arguments
    parser.add_argument(
        '-fq', '--fastq', type=str, help='Correctly demultiplexed/processed FASTQ file containing TCR reads', required=True)
    parser.add_argument(
        '-c', '--chain', type=str, help='TCR chain (a/b/g/d)', required=False)
    parser.add_argument(
        '-br','--bc_read',type=str, help='Which read has bar code (R1,R2)',required=True)
    parser.add_argument(
        '-s', '--suppresssummary', action='store_true', help='Suppress the production of summary data log file', required=False)
    parser.add_argument(
        '-dz', '--dontgzip', action='store_true', help='Stop the output FASTQ files automatically being compressed with gzip', required=False)
    parser.add_argument(
        '-dk', '--dontcheck', action='store_true', help='Skip the FASTQ check', required=False, default=False)  
    parser.add_argument(
        '-dc', '--dontcount', action='store_true', help='Stop Decombinator printing a running count', required=False)
    parser.add_argument(
        '-ex', '--extension', type=str, help='Specify the file extension of the output DCR file. Default = \"n12\"', required=False, default="n12")
    parser.add_argument(
        '-pf', '--prefix', type=str, help='Specify the prefix of the output DCR file. Default = \"dcr_\"', required=False, default="dcr_")
    parser.add_argument(
        '-or', '--orientation', type=str, help='Specify the orientation to search in (forward/reverse/both). Default = reverse', required=False, default="reverse")  
    parser.add_argument(
        '-tg', '--tags', type=str, help='Specify which Decombinator tag set to use (extended or original). Default = extended', required=False, default="extended")
    parser.add_argument(
        '-sp', '--species', type=str, help='Specify which species TCR repertoire the data consists of (human or mouse). Default = human', required=False, default="human")
    parser.add_argument(
        '-N', '--allowNs', action='store_true', help='Whether to allow VJ rearrangements containing ambiguous base calls (\'N\'). Default = False', required=False)
    parser.add_argument(
        '-ln', '--lenthreshold', type=int, help='Acceptable threshold for inter-tag (V to J) sequence length. Default = 130', required=False, default=130)
    parser.add_argument(
        '-tfdir', '--tagfastadir', type=str, help='Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis. \
        Default = \"Decombinator-Tags-FASTAs\".', required=False, default="Decombinator-Tags-FASTAs")
    parser.add_argument(
        '-nbc', '--nobarcoding', action='store_true', help='Option to run Decombinator without barcoding, i.e. so as to run on data produced by any protocol.', required=False)
    parser.add_argument(
        '-bl', '--bclength', type=int, help='Length of barcode sequence, if applicable. Default is set to 42 bp.', required=False, default=42)

    # Collapsinator arguments
    parser.add_argument(
        '-mq', '--minbcQ', type=int, help='Minimum quality score that barcode nucleotides should be to for that rearrangement to be retained. Default = 20.', \
        required=False, default=20)
    parser.add_argument(
        '-bm', '--bcQbelowmin', type=int, help='Number of nucleotides per barcode whose quality score are allowed to be below -mq and still be retained. Default = 1.', \
        required=False, default=1)
    parser.add_argument(
        '-aq', '--avgQthreshold', type=int, help='Average quality threshold that barcode sequences must remain above for rearrangements to be retained. Default = 30', \
        required=False, default=30)
    parser.add_argument(
        '-lv', '--percentlevdist', type=int, help='Percentage Levenshtein distance that is allowed to estimate whether two sequences within a barcode are \
        derived from the same originator molecule. Default = 10', required=False, default=10)  
    parser.add_argument(
        '-bc', '--bcthreshold', type=int, help='Number of sequence edits that are allowed to consider two barcodes to be derived from same originator \
        during clustering. Default = 2.', required=False, default=2)
    parser.add_argument(
        '-di', '--dontcheckinput', action='store_true', help='Override the inputfile sanity check', required=False)
    parser.add_argument(
        '-bd', '--barcodeduplication', action='store_true', help='Optionally output a file containing the final list of clustered barcodes, and their frequencies',\
        required=False)
    parser.add_argument(
        '-pb', '--positionalbarcodes', action='store_true', help='Instead of inferring random barcode sequences from their context relative to spacer sequences, just take the sequence at the default positions. Useful to salvage runs when R2 quality is terrible.',\
        required=False)
    parser.add_argument(
        '-ol', '--oligo', type=str, help='Choose experimental oligo for correct identification of spacers ["M13", "I8","I8_single] (default: M13)',\
        required=True, default="m13")
    parser.add_argument(
        '-wc', '--writeclusters', action='store_true', help='Write cluster data to separate cluster files',\
        required=False, default=False)
    parser.add_argument(
        '-uh', '--UMIhistogram', action='store_true', help='Creates histogram of average UMI cluster sizes',\
        required=False, default=False)

    # CDR3translator arguments
    parser.add_argument('-npf', '--nonproductivefilter', action='store_true', required=False,
                        help='Filter out non-productive reads from the output')
 
    return vars(parser.parse_args())

def create_args_dict(
    fastq: str,
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
    oligo: str = "m13",
    writeclusters: bool = False,
    UMIhistogram: bool = False,
    nonproductivefilter: bool = False
) -> dict:
    
    """
    Creates a function argument dictionary to be used in Decombinator, 
    Collapsinator, and CDR3translator.
    """
    
    return {
        "fastq": fastq,
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
        "nonproductivefilter": nonproductivefilter
    }

# Example usage:
args_dict = create_args_dict('-fq', 'your_fastq_file.fastq', '-br', 'R1')
print(args_dict)


def sort_permissions(fl):
    """
    Need to ensure proper file permissions on output data.
    If users are running pipeline through Docker might otherwise require root access
    :param fl: The file to sort permissions on
    :return: Nothing: script edits permissions where appropriate, if possible
    """

    if oct(os.stat(fl).st_mode)[4:] != '666':
        os.chmod(fl, 0o666)

def write_out_intermediate(data: list, inputargs: dict, suffix: str):
    chain = inputargs["chain"]
    chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}
    filename_id = os.path.basename(inputargs['fastq']).split(".")[0]
    outfilename = f"dcr_{filename_id}" + f"_{chainnams[chain]}" + suffix
    with open(outfilename, 'w') as outfile:
        for line in data:
            outfile.write(", ".join(map(str, line)) + "\n")

    if not inputargs['dontgzip']:
        print("Compressing intermediate output file to", outfilename + ".gz")

        with open(outfilename) as infile, gzip.open(outfilename + '.gz', 'wt') as outfile:
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
    filename_id = os.path.basename(inputargs['fastq']).split(".")[0]
    outfilename = f"dcr_{filename_id}" + f"_{chainnams[chain]}" + suffix
    data.to_csv(f"{outfilename}", sep="\t", index=False)

    if not inputargs['dontgzip']:
        print("Compressing pipeline output file to", outfilename + ".gz")

        with open(outfilename) as infile, gzip.open(outfilename + '.gz', 'wt') as outfile:
            outfile.writelines(infile)
        os.unlink(outfilename)

        outfilenam = outfilename + ".gz"

    else:
        outfilenam = outfilename

    sort_permissions(outfilenam)