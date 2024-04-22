from decombinator import pipeline, io
import pytest
import pathlib
import os
from Bio import BiopythonWarning

@pytest.fixture
def output_dir(tmp_path: pathlib.Path) -> pathlib.Path:
    # Create a temporary output directory
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir

@pytest.fixture
def resource_location() -> pathlib.Path:
    return pathlib.Path("tests/resources")

@pytest.fixture
def reference_file(resource_location: pathlib.Path) -> str:
    # Read reference file from resources folder
    reference_file_path = resource_location / "reference_file.txt"
    with open(reference_file_path, "r") as f:
        reference_content = f.read()
    return reference_content

@pytest.fixture(params=['a', 'b'])
def chain_type(request):
    return request.param

# TODO: Edit decombinator to handle Biopython warnings
@pytest.mark.filterwarnings("ignore::Bio.BiopythonWarning")
def test_race_pipeline(
        output_dir: pathlib.Path,
        resource_location: pathlib.Path,
        chain_type: str
    ) -> None:

    filename: str = "TINY_1.fq.gz"
    args = io.create_args_dict(
        fastq=str((
            resource_location / filename
            ).resolve()
            )
        ,
        chain=chain_type,
        bc_read="R2",
        dontgzip=True,
        outpath=f"{output_dir}{os.sep}",
    )
    pipeline.run(args)