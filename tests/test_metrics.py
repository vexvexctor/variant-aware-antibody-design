"""Sign convention and cluster-bootstrap unit tests."""
import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from analysis.metrics import (beat_native_rate, beats_native, best_of_pool,
                              paired_targets, worst_heldout_margin)
from analysis.bootstrap import cluster_bootstrap, paired_delta_bootstrap


class TestMargins(unittest.TestCase):
    def test_worst_is_the_maximum_margin(self):
        self.assertEqual(worst_heldout_margin([-2.0, -0.5, -3.0]), -0.5)

    def test_beats_native_requires_every_variant_below_zero(self):
        self.assertTrue(beats_native([-2.0, -0.5, -3.0]))
        self.assertFalse(beats_native([-2.0, 0.1, -3.0]))   # one bad variant is enough

    def test_best_of_pool_picks_least_bad_worst_case(self):
        self.assertEqual(best_of_pool([[-0.1, 2.0], [-1.0, -0.5], [0.0, 0.0]]), -0.5)

    def test_rate(self):
        self.assertAlmostEqual(beat_native_rate({"a": -1.0, "b": 0.5, "c": -0.2}), 2 / 3)

    def test_paired_intersects(self):
        self.assertEqual(paired_targets({"a": 1, "b": 2}, {"b": 3, "c": 4}), ["b"])


class TestBootstrap(unittest.TestCase):
    def test_is_seeded_and_reproducible(self):
        v, c = [1.0, 2.0, 3.0, 4.0], ["x", "x", "y", "y"]
        self.assertEqual(cluster_bootstrap(v, c, seed=0), cluster_bootstrap(v, c, seed=0))

    def test_resamples_clusters_not_rows(self):
        """Same 20 rows: 2 coarse clusters must give a WIDER interval than 20 singletons,
        because the effective sample size is the number of clusters, not the row count."""
        v = [0.0] * 10 + [1.0] * 10                       # clusters must differ in mean
        _, lo_few, hi_few = cluster_bootstrap(v, ["a"] * 10 + ["b"] * 10, seed=1)
        _, lo_many, hi_many = cluster_bootstrap(v, [str(i) for i in range(20)], seed=1)
        self.assertGreater(hi_few - lo_few, hi_many - lo_many)

    def test_paired_delta_sign(self):
        a = {"t1": 0.5, "t2": 0.5}        # neither beats native
        b = {"t1": -0.5, "t2": -0.5}      # both do
        d, _, _ = paired_delta_bootstrap(a, b, {"t1": "c1", "t2": "c2"}, draws=200)
        self.assertAlmostEqual(d, 100.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
