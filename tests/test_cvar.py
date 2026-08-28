"""CVaR sign/tail tests. A sign error here would materially change the method."""
import math, os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from steering.cvar import aggregate, cvar_margin, cvar_reward, objective, tail_size


class TestTailSize(unittest.TestCase):
    def test_matches_ceil_rule(self):
        self.assertEqual(tail_size(6, 0.2), 2)      # ceil(1.2)
        self.assertEqual(tail_size(30, 0.2), 6)
        self.assertEqual(tail_size(6, 1.0), 6)
    def test_never_empty(self):
        self.assertEqual(tail_size(6, 0.0), 1)
        self.assertEqual(tail_size(1, 0.2), 1)


class TestRewardSpace(unittest.TestCase):
    R = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    def test_cvar_takes_the_lowest_rewards(self):
        self.assertAlmostEqual(cvar_reward(self.R, 0.2), -1.5)   # mean(-2,-1)
    def test_alpha_one_is_the_mean(self):
        self.assertAlmostEqual(aggregate(self.R, "cvar", 1.0), aggregate(self.R, "mean"))
    def test_small_alpha_is_the_worst_case(self):
        self.assertAlmostEqual(aggregate(self.R, "cvar", 1e-9), aggregate(self.R, "worst"))
    def test_order_invariant(self):
        self.assertAlmostEqual(cvar_reward(self.R, .2), cvar_reward(list(reversed(self.R)), .2))
    def test_monotone_in_alpha(self):
        vals = [cvar_reward(self.R, a) for a in (0.1, 0.2, 0.5, 1.0)]
        self.assertEqual(vals, sorted(vals))       # more lenient tail => higher reward


class TestMarginSpace(unittest.TestCase):
    M = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
    def test_cvar_margin_takes_the_largest_margins(self):
        # larger margin = worse, so the bad tail is the UPPER tail
        self.assertAlmostEqual(cvar_margin(self.M, 0.2), 1.5)    # mean(1,2)
    def test_objective_is_negated_margin_cvar(self):
        self.assertAlmostEqual(objective(self.M, 0.2), -cvar_margin(self.M, 0.2))
    def test_reward_and_margin_views_agree(self):
        rewards = [-m for m in self.M]
        self.assertAlmostEqual(objective(self.M, 0.2), cvar_reward(rewards, 0.2))
    def test_negative_margin_means_better_than_native(self):
        # a design beating native on every variant must score a positive objective
        self.assertGreater(objective([-1.0] * 6, 0.2), 0.0)


class TestEquivalenceWithSteeringLoop(unittest.TestCase):
    """The in-loop copy in variant_aware.py must agree with this reference."""
    def test_matches_inline_implementation(self):
        try:
            import torch
        except ImportError:
            self.skipTest("torch not installed")
        for vals in ([-2., -1., 0., 1., 2., 3.], [0.5] * 6, [-5., 4., 4., 4., 4., 4.]):
            for alpha in (0.1, 0.2, 0.5, 1.0):
                v = torch.tensor(vals, dtype=torch.float)
                k = max(1, math.ceil(alpha * v.numel()))
                inline = float(torch.sort(v).values[:k].mean())   # verbatim from variant_aware.py
                self.assertAlmostEqual(inline, aggregate(vals, "cvar", alpha), places=5)


if __name__ == "__main__":
    unittest.main(verbosity=2)
