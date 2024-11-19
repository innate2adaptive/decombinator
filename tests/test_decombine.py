import pathlib
import time
from typing import Any

import pytest

from decombinator import decombine


class TestEmptyFq:
    """
    Testing that empty FASTQ is handled gracefully with logging.
    """

    @pytest.fixture(scope="class")
    def output_dir(
        self, tmp_path_factory: pytest.TempPathFactory
    ) -> pathlib.Path:
        output_dir = tmp_path_factory.mktemp("output")
        return output_dir

    @pytest.fixture
    def empty_filepath(self, output_dir: pathlib.Path) -> pathlib.Path:
        return output_dir / "empty_merge.fq"

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
            "outpath": str(output_dir.resolve()) + "/",
            "suppresssummary": False,
            "tags": "extended",
            "species": "human",
            "tagfastadir": "tests/resources/Decombinator-Tags-FASTAs",
        }

    @pytest.fixture
    def test_empty_fq(
        self, empty_file: None, pipe_args: dict[str, Any]
    ) -> None:
        with pytest.raises(ValueError):
            decombine.decombinator(pipe_args)

    @pytest.fixture
    def expected_log(self) -> str:
        return "OutputFile," + "empty_alpha" + "\nNumberReadsInput," + "0\n"

    @pytest.fixture
    def test_empty_log(
        self, output_dir: pathlib.Path, test_empty_fq: None, expected_log: str
    ) -> None:
        date = time.strftime("%Y_%m_%d")
        logfile = (
            output_dir
            / "Logs"
            / f"{date}_alpha_empty_merge_Decombinator_Summary.csv"
        )
        with logfile.open() as log:
            assert log.read() == expected_log

    @pytest.fixture
    def test_empty_fq_repeat(
        self, pipe_args: dict[str, Any], test_empty_log: None
    ) -> None:
        with pytest.raises(ValueError):
            decombine.decombinator(pipe_args)

    # TODO: Test this logic independently. Due to logging design in Decombinator
    # at present this must be tested in the empty fq case specifically
    def test_2nd_log(
        self,
        output_dir: pathlib.Path,
        test_empty_fq: None,
        expected_log: str,
        test_empty_fq_repeat: None,
    ) -> None:
        date = time.strftime("%Y_%m_%d")
        logfile = (
            output_dir
            / "Logs"
            / f"{date}_alpha_empty_merge_Decombinator_Summary2.csv"
        )
        with logfile.open() as log:
            assert log.read() == expected_log
