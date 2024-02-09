from Decombinator import decombinator
from Collapsinator import collapsinator
from CDR3translator import cdr3translator
from dcr_utilities import args, write_out_translated, write_out_intermediate

from datetime import datetime
startTime = datetime.now()

if __name__ == '__main__':

    inputargs = args()

    # Run pipline, ovewriting data after each function call to save memory
    data = decombinator(inputargs)
    write_out_intermediate(data, inputargs, ".n12")
    print("Decombinator complete...")

    data = collapsinator(data, inputargs)
    write_out_intermediate(data, inputargs, ".freq")
    print("Collapsinator complete...")

    data = cdr3translator(data, inputargs)
    print("CDR3translator complete...")

    write_out_translated(data, inputargs)
    print(f"Pipeline complete in {datetime.now() - startTime}")