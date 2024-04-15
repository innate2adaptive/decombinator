from decombinator import pipeline, io
import pytest
import pathlib

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

def test_pipeline_integration_race(
        output_dir: pathlib.Path,
        resource_location: pathlib.Path
    ) -> None:

    # TODO use a fixture to perform this test for both A and B
    filename: str = "TINY_1.fq.gz"
    args = io.create_args_dict(
        fastq=(
        resource_location / filename
        ).resolve().as_posix(),
        chain="B",
        bc_read="R2",
        dontgzip=True,
        prefix="output/dcr_",
    )
    pipeline.run(args)

