import collections as coll
import pathlib
import typing

import pytest

from decombinator import collapse


class TestClusterUMIs:
    def test_no_umis(self):
        with pytest.raises(ValueError):
            collapse.cluster_UMIs(coll.defaultdict(list), {}, 0, 0, False)

    @pytest.fixture
    def barcode_dcretc_list(self):
        return {
            "AAAA|0|AAAA": ["AAAA"],
            "GGGG|0|GGGG": ["GGGG"],
            "AAAG|0|AAAG": ["AAAG"],
            "AAAA|1|GGGG": ["GGGG"],
        }

    def test_create_cluster_structures(self, barcode_dcretc_list):
        num_initial_groups, barcode_dcretc_list, umi_protoseq_tuple = (
            collapse.create_clustering_objs(barcode_dcretc_list)
        )

        assert num_initial_groups == 4
        assert barcode_dcretc_list == [
            ("AAAA|0|AAAA", ["AAAA"]),
            ("GGGG|0|GGGG", ["GGGG"]),
            ("AAAG|0|AAAG", ["AAAG"]),
            ("AAAA|1|GGGG", ["GGGG"]),
        ]
        assert umi_protoseq_tuple == [
            ("AAAA", "AAAA"),
            ("GGGG", "GGGG"),
            ("AAAG", "AAAG"),
            ("AAAA", "GGGG"),
        ]

    def test_merge_order(self, barcode_dcretc_list):
        clusters = collapse.cluster_UMIs(
            barcode_dcretc_list, {"writeclusters": False}, 2, 0.25, True
        )

        assert clusters == {
            "AAAA|0|AAAA": ["AAAA", "AAAG"],
            "GGGG|0|GGGG": ["GGGG"],
            "AAAA|1|GGGG": ["GGGG"],
        }


class TestGetBarcodePositions:

    @pytest.fixture
    def counter(self):
        return coll.Counter()

    def test_m13(self, counter):
        m13 = "GTCGTGACTGGGAAAACCCTGG"
        i8 = "GTCGTGAT"
        bcseq = "GTCGTGACTGGGAAAACCCTGGTTTCCGGTCGTGATAAAGTG"
        inputargs = {
            "oligo": "m13",
            "allowNs": False,
        }

        assert collapse.get_barcode_positions(bcseq, inputargs, counter) == [
            len(m13),
            len(m13) + 6,
            len(m13) + 6 + len(i8),
            len(m13) + 6 + len(i8) + 6,
        ]

    def test_i8(self, counter):
        i8 = "GTCGTGAT"
        bcseq = "GTCGTGATTTTCCGGTCGTGATAAAGTG"
        inputargs = {
            "oligo": "i8",
            "allowNs": False,
        }

        assert collapse.get_barcode_positions(bcseq, inputargs, counter) == [
            len(i8),
            len(i8) + 6,
            len(i8) + 6 + len(i8),
            len(i8) + 6 + len(i8) + 6,
        ]

    def test_i8_single(self, counter):
        i8 = "ATCACGAC"
        bcseq = "GAAGCTATCACGACATCACTAC"
        inputargs = {
            "oligo": "i8_single",
            "allowNs": False,
        }

        assert collapse.get_barcode_positions(bcseq, inputargs, counter) == [
            0,
            6,
            6 + len(i8),
            6 + len(i8) + 6,
        ]

    def test_nebio(self, counter):
        bcseq = "CGGGCTTGGTATCGGCCGATCTACGGG"
        inputargs = {
            "oligo": "nebio",
            "allowNs": False,
        }

        assert collapse.get_barcode_positions(bcseq, inputargs, counter) == [
            0,
            17,
        ]

    def test_takara(self, counter):
        bcseq = "CTCGTTAGGTTCGTACGGGGATTGCA"
        inputargs = {
            "oligo": "takara",
            "allowNs": False,
        }

        assert collapse.get_barcode_positions(bcseq, inputargs, counter) == [
            0,
            12,
        ]


class TestFindFirstSpacer:

    def test_m13(self):
        oligo = {
            "spcr1": "GTCGTGACTGGGAAAACCCTGG",
        }
        seq = "GTCGTGACTGGGAAAACCCTGGTTTCCGGTCGTGATAAAGTG"
        oligo_start = 0
        allowance = 10
        oligo_end = allowance + len(oligo["spcr1"])

        assert collapse.findFirstSpacer(
            oligo, seq, oligo_start, oligo_end
        ) == [oligo["spcr1"]]

    def test_i8(self):
        oligo = {
            "spcr1": "GTCGTGAT",
        }
        seq = "GTCGTGATTTTCCGGTCGTGATAAAGTG"
        oligo_start = 0
        allowance = 10
        oligo_end = allowance + len(oligo["spcr1"])

        assert collapse.findFirstSpacer(
            oligo, seq, oligo_start, oligo_end
        ) == [oligo["spcr1"]]

    def test_i8_single(self):
        oligo = {
            "spcr1": "ATCACGAC",
        }
        seq = "GAAGCTATCACGACATCACTAC"
        oligo_start = 0
        allowance = 10
        oligo_end = allowance + len(oligo["spcr1"])

        assert collapse.findFirstSpacer(
            oligo, seq, oligo_start, oligo_end
        ) == [oligo["spcr1"]]

    def test_nebio(self):
        oligo = {
            "spcr1": "TACGGG",
        }
        seq = "CGGGCTTGGTATCGGCCGATCTACGGG"
        oligo_start = 18
        oligo_end = oligo_start + 10

        assert collapse.findFirstSpacer(
            oligo, seq, oligo_start, oligo_end
        ) == [oligo["spcr1"]]

    def test_takara(self):
        oligo = {
            "spcr1": "GTACGGG",
        }
        seq = "CTCGTTAGGTTCGTACGGGGATTGCA"
        oligo_start = 0
        oligo_end = 19

        assert collapse.findFirstSpacer(
            oligo, seq, oligo_start, oligo_end
        ) == [oligo["spcr1"]]


class TestReadInData:

    collapse.counts = coll.Counter()

    @pytest.fixture
    def blank_input(self) -> list[None]:
        return []

    @pytest.fixture
    def pipe_args(self) -> dict[str, typing.Union[str, int]]:
        return {
            "command": "pipeline",
            "lenthreshold": 130,
            "minbcQ": 20,
            "bcQbelowmin": 1,
            "avgQthreshold": 30,
            "oligo": "M13",
            "sampling_analysis": False,
        }

    def test_no_dcr(self, blank_input, pipe_args):
        with pytest.raises(ValueError):
            collapse.read_in_data(
                blank_input, pipe_args, None, None, None, None
            )

    @pytest.fixture
    def valid_input(self) -> list[list[str]]:
        return [
            [
                "15",
                "4",
                "1",
                "7",
                "CCCCCAGGGGGCTC",
                "LH00409:259:22JJCFLT4:8:1149:38410:17375",
                "ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGCCCCCAGGGGGCTCCAGCCCCAGCATTTTGGTGATGGGACTCGACTC",
                "IIIIIIIIIIIIIII-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                "GTCGTGACTGGGAAAACCCTGGCACCCGGTCGTGATCTGACT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                "GACACAACTCTCCCCAGAGAAGGTGGTGTGA",
            ],
            [  # Final base of junction is changed, should still collapse
                "15",
                "4",
                "1",
                "7",
                "CCCCCAGGGGGCTG",
                "LH00409:259:22JJCFLT4:8:1149:38410:17375",
                "ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGCCCCCAGGGGGCTGCAGCCCCAGCATTTTGGTGATGGGACTCGACTC",
                "IIIIIIIIIIIIIII-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                "GTCGTGACTGGGAAAACCCTGGCACCCGGTCGTGATCTGACT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                "GACACAACTCTCCCCAGAGAAGGTGGTGTGA",
            ],
            [  # Same TCR but with different junction region, making TCR 20% different
                "15",
                "4",
                "1",
                "7",
                "AAAAAAAAAAAAAA",
                "LH00409:259:22JJCFLT4:8:1149:38410:17375",
                "ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGAAAAAAAAAAAAAACAGCCCCAGCATTTTGGTGATGGGACTCGACTC",
                "IIIIIIIIIIIIIII-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                "GTCGTGACTGGGAAAACCCTGGCACCCGGTCGTGATCTGACT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                "GACACAACTCTCCCCAGAGAAGGTGGTGTGA",
            ],
        ]

    @pytest.fixture
    def barcode_quality_parameters(
        self, pipe_args: dict[str, int]
    ) -> list[int]:
        return [
            pipe_args["minbcQ"],
            pipe_args["bcQbelowmin"],
            pipe_args["avgQthreshold"],
        ]

    def test_barcode_collision(
        self,
        valid_input: list[list[str]],
        pipe_args: dict[str, str | int],
        barcode_quality_parameters: list[int],
    ):
        barcode_dcretc = collapse.read_in_data(
            valid_input,
            pipe_args,
            barcode_quality_parameters,
            lev_threshold_fraction=0.1,
            dont_count=False,
            opener=open,
        )
        assert barcode_dcretc == {
            "CACCCGCTGACT|0|ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGCCCCCAGGGGGCTCCAGCCCCAGCATTTTGGTGATGGGACTCGACTC": [
                "['15', '4', '1', '7', 'CCCCCAGGGGGCTC']|ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGCCCCCAGGGGGCTCCAGCCCCAGCATTTTGGTGATGGGACTCGACTC|IIIIIIIIIIIIIII-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII|LH00409:259:22JJCFLT4:8:1149:38410:17375",
                "['15', '4', '1', '7', 'CCCCCAGGGGGCTG']|ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGCCCCCAGGGGGCTGCAGCCCCAGCATTTTGGTGATGGGACTCGACTC|IIIIIIIIIIIIIII-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII|LH00409:259:22JJCFLT4:8:1149:38410:17375",
            ],
            "CACCCGCTGACT|1|ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGAAAAAAAAAAAAAACAGCCCCAGCATTTTGGTGATGGGACTCGACTC": [
                "['15', '4', '1', '7', 'AAAAAAAAAAAAAA']|ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGAAAAAAAAAAAAAACAGCCCCAGCATTTTGGTGATGGGACTCGACTC|IIIIIIIIIIIIIII-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII|LH00409:259:22JJCFLT4:8:1149:38410:17375"
            ],
        }

    def test_barcode_collision_no_tcr_check(
        self,
        valid_input: list[list[str]],
        pipe_args: dict[str, str | int],
        barcode_quality_parameters: list[int],
    ):
        barcode_dcretc = collapse.read_in_data(
            valid_input,
            pipe_args,
            barcode_quality_parameters,
            lev_threshold_fraction=1,
            dont_count=False,
            opener=open,
        )
        assert barcode_dcretc == {
            "CACCCGCTGACT|0|ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGCCCCCAGGGGGCTCCAGCCCCAGCATTTTGGTGATGGGACTCGACTC": [
                "['15', '4', '1', '7', 'CCCCCAGGGGGCTC']|ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGCCCCCAGGGGGCTCCAGCCCCAGCATTTTGGTGATGGGACTCGACTC|IIIIIIIIIIIIIII-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII|LH00409:259:22JJCFLT4:8:1149:38410:17375",
                "['15', '4', '1', '7', 'CCCCCAGGGGGCTG']|ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGCCCCCAGGGGGCTGCAGCCCCAGCATTTTGGTGATGGGACTCGACTC|IIIIIIIIIIIIIII-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII|LH00409:259:22JJCFLT4:8:1149:38410:17375",
                "['15', '4', '1', '7', 'AAAAAAAAAAAAAA']|ATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGAAAAAAAAAAAAAACAGCCCCAGCATTTTGGTGATGGGACTCGACTC|IIIIIIIIIIIIIII-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIII-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII|LH00409:259:22JJCFLT4:8:1149:38410:17375",
            ],
        }


class TestCheckDcrFile:

    @pytest.fixture(scope="class")
    def output_dir(
        self, tmp_path_factory: pytest.TempPathFactory
    ) -> pathlib.Path:
        output_dir = tmp_path_factory.mktemp("output")
        return output_dir

    @pytest.fixture
    def empty_filepath(self, output_dir: pathlib.Path) -> pathlib.Path:
        return output_dir / "empty.n12"

    @pytest.fixture
    def empty_file(self, empty_filepath: pathlib.Path) -> None:
        output = ""
        empty_filepath.write_text(output)

    def test_empty_n12(self, empty_file, empty_filepath: pathlib.Path) -> None:
        with pytest.raises(ValueError):
            collapse.check_dcr_file(empty_filepath, open)
