from decombinator import collapse
import collections as coll
import pytest


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
            barcode_dcretc_list, {"writeclusters": False}, 2, 25, True
        )

        assert clusters == {
            "AAAA|0|AAAA": ["AAAA", "AAAG"],
            "GGGG|0|GGGG": ["GGGG"],
            "AAAA|1|GGGG": ["GGGG"],
        }
