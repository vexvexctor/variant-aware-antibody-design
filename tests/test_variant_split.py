"""The design/held-out split must be disjoint -- the benchmark's whole validity."""
import csv, os, unittest
ROOT = os.path.join(os.path.dirname(__file__), "..")
M = os.path.join(ROOT, "data", "manifests")
S = os.path.join(ROOT, "data", "splits")


def _by_target(path, key="variant_name"):
    out = {}
    for r in csv.DictReader(open(path)):
        out.setdefault(r["target"], set()).add(r[key])
    return out


class TestSplit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.design = _by_target(os.path.join(S, "design_variants.csv"))
        cls.heldout = _by_target(os.path.join(S, "heldout_variants.csv"))

    def test_no_overlap_on_any_target(self):
        for t, d in self.design.items():
            self.assertEqual(d & self.heldout.get(t, set()), set(),
                             f"design/held-out overlap on {t}")

    def test_panel_sizes(self):
        self.assertTrue(all(len(v) == 6 for v in self.design.values()))
        self.assertTrue(all(len(v) == 24 for v in self.heldout.values()))

    def test_target_count(self):
        self.assertEqual(len(self.design), 200)
        self.assertEqual(len(self.heldout), 200)

    def test_every_target_has_both_panels(self):
        self.assertEqual(set(self.design), set(self.heldout))


class TestClusters(unittest.TestCase):
    def setUp(self):
        self.rows = list(csv.DictReader(open(os.path.join(M, "antigen_clusters.csv"))))

    def test_benchmark_labels_are_the_fixed_172(self):
        self.assertEqual(len({r["benchmark_cluster"] for r in self.rows}), 172)

    def test_full_aacdb_labels_are_the_171(self):
        """171 != 172 by construction: different clustering scopes. See data/README.md."""
        self.assertEqual(len({r["full_aacdb_cluster"] for r in self.rows if r["full_aacdb_cluster"]}), 171)

    def test_every_target_has_a_benchmark_cluster(self):
        self.assertTrue(all(r["benchmark_cluster"] for r in self.rows))


if __name__ == "__main__":
    unittest.main(verbosity=2)
