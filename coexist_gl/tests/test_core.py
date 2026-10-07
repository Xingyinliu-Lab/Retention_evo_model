from __future__ import annotations

import itertools
import math
import unittest
from io import StringIO

import numpy as np
from Bio import Phylo
from scipy.linalg import expm

from coexist_gl.clustering import clustering_test
from coexist_gl.abstract_model import AbstractTreeLikelihood
from coexist_gl.model import StateTreeLikelihood, root_prior, transition_matrix
from coexist_gl.species_cuts import discover_topology_cut


class CoreTests(unittest.TestCase):
    def test_r1_transition_identity(self):
        for A in (0.01, 0.5, 2.0, 20.0):
            expected = expm(np.array([[-A, A], [2.0, -2.0]]))
            observed = transition_matrix(A, 1.0)
            np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(observed.sum(axis=1), 1.0, atol=1e-12)

    def test_pruning_matches_direct_enumeration(self):
        tree = Phylo.read(StringIO("((a,b),c);"), "newick")
        states = {"a": 1, "b": 1, "c": 0}
        model = StateTreeLikelihood(tree, states)
        A, r = 0.7, 3.0
        observed = math.exp(model.evaluate(A, r).log_likelihood)
        P = transition_matrix(A, r)
        pi = root_prior(A)
        direct = 0.0
        for root_state, internal_state in itertools.product(range(2), repeat=2):
            direct += (
                pi[root_state]
                * P[root_state, internal_state]
                * P[root_state, states["c"]]
                * P[internal_state, states["a"]]
                * P[internal_state, states["b"]]
            )
        self.assertAlmostEqual(observed, direct, places=12)

    def test_edge_expectations_sum_to_edge_count(self):
        tree = Phylo.read(StringIO("((a,b),(c,d));"), "newick")
        model = StateTreeLikelihood(tree, {"a": 0, "b": 1, "c": 0, "d": 1})
        result = model.evaluate(0.8, 2.0, with_expectations=True)
        self.assertAlmostEqual(result.expected["expected_total_edges"], 6.0, places=10)

    def test_clustering_output(self):
        tree = Phylo.read(StringIO("(((a,b),(c,d)),((e,f),(g,h)));"), "newick")
        states = {name: int(name in {"a", "b", "c", "d"}) for name in "abcdefgh"}
        result = clustering_test(tree, states, permutations=999, seed=7)
        self.assertEqual(result["fitch_steps_observed"], 1)
        self.assertLess(result["p_cluster"], 0.05)

    def test_species_cut_is_state_blind_and_mass_preserving(self):
        tree = Phylo.read(StringIO("(((a,b),(c,d)),((e,f),(g,h)));"), "newick")
        states_a = {name: int(name in {"a", "b", "c"}) for name in "abcdefgh"}
        states_b = {name: 1 - states_a[name] for name in "abcdefgh"}
        cut_a = discover_topology_cut(tree, states_a, 4)
        cut_b = discover_topology_cut(tree, states_b, 4)
        membership_a = [(x["species_block"], x["host_id"]) for x in cut_a.members]
        membership_b = [(x["species_block"], x["host_id"]) for x in cut_b.members]
        self.assertEqual(membership_a, membership_b)
        self.assertEqual(len(membership_a), 8)
        self.assertTrue(math.isfinite(AbstractTreeLikelihood(cut_a.root).evaluate(0.5, 2.0, 10.0).log_likelihood))


if __name__ == "__main__":
    unittest.main()
