from Decombinator import decombinator
from Collapsinator import collapsinator
from CDR3translator import cdr3translator
from dcr_utilities import args, write_out

from datetime import datetime
startTime = datetime.now()

if __name__ == '__main__':

    inputargs = vars(args())

    # Run pipline, ovewriting data after each function call to save memory
    data = decombinator(inputargs)
    print("Decombinator complete..")

    data = collapsinator(data, inputargs)
    print("Collapsinator complete...")

    data = cdr3translator(data, inputargs)
    print("CDR3translator complete...")

    write_out(data, inputargs)
    print(f"Pipeline complete in {datetime.now() - startTime}")