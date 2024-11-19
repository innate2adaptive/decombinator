import pathlib
from typing import Any

import pytest

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
    def pipe_args(
        self, output_dir: pathlib.Path, empty_filepath: pathlib.Path
    ) -> dict[str, Any]:
        return {
            "command": "pipeline",
            "infile": str(empty_filepath.resolve()),
            "dontcheck": False,
            "chain": "a",
            "outpath": str(output_dir.resolve()),
        }

    @pytest.fixture
    def test_empty_fq(
        self, empty_file: None, pipe_args: dict[str, Any]
    ) -> None:
        with pytest.raises(ValueError):
            decombine.decombinator(pipe_args)

    def test_empty_log(
        self, output_dir: pathlib.Path, test_empty_fq: None
    ) -> None:
        logfile = output_dir / "Logs" / "empty_Decombinator_Summary.csv"
        with logfile.open() as log:
            inputs = log.readline(21)
            assert inputs == "NumberReadsInput,0"
