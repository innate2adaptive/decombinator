from decombinator import collapse
import collections as coll
import pytest


class TestClusterUMIs:
    def test_no_umis(self):
        with pytest.raises(ValueError):
            collapse.cluster_UMIs(coll.defaultdict(list), {}, 0, 0, False)
