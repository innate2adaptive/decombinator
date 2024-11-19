import pytest
import pathlib
from decombinator import decombine

class TestEmptyFq:

    @pytest.fixture(scope="class")
    def output_dir(
        self, tmp_path_factory: pytest.TempPathFactory
    ) -> pathlib.Path:
        output_dir = tmp_path_factory.mktemp("output")
        return output_dir

    @pytest.fixture
    def empty_filepath(self, output_dir: pathlib.Path) -> pathlib.Path:
        return output_dir / "empty.fq"

    @pytest.fixture
    def empty_file(self, empty_filepath: pathlib.Path) -> None:
        output = ""
        empty_filepath.write_text(output)

    @pytest.fixture
    def pipe_args(self, empty_filepath):
        return {
            "command": "pipeline",
            "infile": str(empty_filepath.resolve()),
            "dontcheck": False
        }

    def test_empty_n12(self, empty_file, empty_filepath: pathlib.Path, pipe_args) -> None:
        with pytest.raises(ValueError):
            decombine.decombinator(pipe_args)
 
