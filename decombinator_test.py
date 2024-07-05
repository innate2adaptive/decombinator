from src.decombinator import decombine, translate, io

if __name__ == '__main__':
    inputargs = io.cli_args()
    data = decombine.decombinator(inputargs)
    data = translate.cdr3translator(data, inputargs)
    io.write_out_translated(data, inputargs)