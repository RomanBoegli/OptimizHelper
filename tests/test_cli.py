"""
CLI integration tests.  Each test invokes a command via Click's CliRunner
and asserts on exit code and key output strings.

Test data lives in tests/data/ (created by the project setup script).
Commands that write output files (plot, branchbound, drawgraph) run inside
an isolated temporary filesystem so they do not pollute the project root.
"""
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from main import main  # conftest.py sets matplotlib backend before this import

# ── paths ─────────────────────────────────────────────────────────────────────
DATA = Path(__file__).parent / 'data'
MATANALYSIS = str(DATA / 'matanalysis.ods')
AB          = str(DATA / 'Ab.ods')
ABC_SIMPLEX = str(DATA / 'Abc_simplex.ods')
ABC_BB      = str(DATA / 'Abc_branchbound.ods')
ADJ_CSV     = str(DATA / 'graph_adjacency.csv')
EDGE_CSV    = str(DATA / 'graph_edgelist.csv')

GRAPHVIZ_AVAILABLE = shutil.which('dot') is not None


@pytest.fixture
def runner():
    return CliRunner()


# ── Part 1a ───────────────────────────────────────────────────────────────────

class TestMatanalysis:
    def test_basic(self, runner):
        r = runner.invoke(main, ['matanalysis', MATANALYSIS])
        assert r.exit_code == 0, r.output
        assert 'rank' in r.output

    def test_pretty(self, runner):
        r = runner.invoke(main, ['matanalysis', '--pretty', MATANALYSIS])
        assert r.exit_code == 0, r.output
        assert 'rank' in r.output


class TestHyperplanes:
    def test_basic(self, runner):
        r = runner.invoke(main, ['hyperplanes', AB])
        assert r.exit_code == 0, r.output
        assert 'xBi' in r.output

    def test_pretty(self, runner):
        r = runner.invoke(main, ['hyperplanes', '--pretty', AB])
        assert r.exit_code == 0, r.output
        assert 'xBi' in r.output


class TestSimplex:
    def test_finds_optimum(self, runner):
        # start from vertex (0,0) = basic sel rows 1,2 (1-indexed)
        r = runner.invoke(main, ['simplex', ABC_SIMPLEX, '1', '2'])
        assert r.exit_code == 0, r.output
        assert 'OPTIMAL' in r.output

    def test_reports_optimal_value(self, runner):
        r = runner.invoke(main, ['simplex', ABC_SIMPLEX, '1', '2'])
        assert r.exit_code == 0, r.output
        assert 'optimal_value' in r.output

    def test_invalid_basic_sel_dimension(self, runner):
        r = runner.invoke(main, ['simplex', ABC_SIMPLEX, '1'])
        assert r.exit_code == 0, r.output
        assert 'does not align' in r.output


class TestPlot:
    def test_creates_png(self, runner):
        with runner.isolated_filesystem():
            r = runner.invoke(main, ['plot', AB])
            assert r.exit_code == 0, r.output
            assert 'plot.png' in r.output

    def test_custom_axis_limits(self, runner):
        with runner.isolated_filesystem():
            r = runner.invoke(main, ['plot', AB, '--xlim', '-5', '5', '--ylim', '-5', '5'])
            assert r.exit_code == 0, r.output


@pytest.mark.skipif(not GRAPHVIZ_AVAILABLE, reason='graphviz (dot) not installed')
class TestBranchBound:
    def test_basic(self, runner):
        with runner.isolated_filesystem():
            r = runner.invoke(main, ['branchbound', ABC_BB])
            assert r.exit_code == 0, r.output
            assert 'tree.png' in r.output

    def test_round_down_first(self, runner):
        with runner.isolated_filesystem():
            r = runner.invoke(main, ['branchbound', ABC_BB, '--rounddownfirst'])
            assert r.exit_code == 0, r.output

    def test_knapsack(self, runner):
        with runner.isolated_filesystem():
            r = runner.invoke(main, ['branchbound', '--knapsack', ABC_BB])
            assert r.exit_code == 0, r.output


# ── Part 2a ───────────────────────────────────────────────────────────────────

class TestDiffbeauty:
    def test_single_var(self, runner):
        r = runner.invoke(main, ['diffbeauty', 'x**2'])
        assert r.exit_code == 0, r.output
        assert '2' in r.output  # d/dx x² = 2x

    def test_wrt_specific(self, runner):
        r = runner.invoke(main, ['diffbeauty', 'x**2 + y**2', '--wrt', 'x'])
        assert r.exit_code == 0, r.output
        assert '2' in r.output

    def test_wrt_all(self, runner):
        r = runner.invoke(main, ['diffbeauty', 'x**2 + y**2', '--wrt', 'all'])
        assert r.exit_code == 0, r.output
        assert r.output.count('2') >= 2  # partial w.r.t. x and y both yield 2·var


class TestDifftree:
    def test_multivar(self, runner):
        r = runner.invoke(main, ['difftree', 'x**2 + y'])
        assert r.exit_code == 0, r.output
        assert 'f' in r.output  # tree nodes are labelled f…

    def test_single_var(self, runner):
        r = runner.invoke(main, ['difftree', 'x**3'])
        assert r.exit_code == 0, r.output


class TestEvaluate:
    def test_two_vars(self, runner):
        # f(x,y) = x²+y, f(2,3) = 4+3 = 7
        r = runner.invoke(main, ['evaluate', 'x**2 + y', '2', '3'])
        assert r.exit_code == 0, r.output
        assert '7' in r.output

    def test_one_var(self, runner):
        r = runner.invoke(main, ['evaluate', 'x**2', '3'])
        assert r.exit_code == 0, r.output
        assert '9' in r.output

    def test_missing_value_error(self, runner):
        r = runner.invoke(main, ['evaluate', 'x**2 + y', '2'])
        assert r.exit_code == 0, r.output
        assert 'Missing' in r.output


class TestGradient:
    def test_basic(self, runner):
        r = runner.invoke(main, ['gradient', 'x**2 + y**2'])
        assert r.exit_code == 0, r.output

    def test_with_substitution(self, runner):
        # ∇f at (1,1): [2, 2]
        r = runner.invoke(main, ['gradient', 'x**2 + y**2', '-s', 'x', '1.0', '-s', 'y', '1.0'])
        assert r.exit_code == 0, r.output
        assert '2' in r.output

    def test_pretty(self, runner):
        r = runner.invoke(main, ['gradient', '--pretty', 'x**2 + y**2'])
        assert r.exit_code == 0, r.output


class TestHessian:
    def test_basic(self, runner):
        r = runner.invoke(main, ['hessian', 'x**2 + y**2'])
        assert r.exit_code == 0, r.output

    def test_determinant(self, runner):
        # H(x²+y²) = diag(2,2), det = 4
        r = runner.invoke(main, ['hessian', 'x**2 + y**2', '--det'])
        assert r.exit_code == 0, r.output
        assert '4' in r.output

    def test_with_substitution(self, runner):
        r = runner.invoke(main, ['hessian', 'x**2 + y**2', '-s', 'x', '1.0', '-s', 'y', '1.0'])
        assert r.exit_code == 0, r.output

    def test_pretty(self, runner):
        r = runner.invoke(main, ['hessian', '--pretty', 'x**2 + y**2'])
        assert r.exit_code == 0, r.output


class TestNewton:
    def test_one_step(self, runner):
        # f=x²+y², start (1,1): H⁻¹=0.5I, ∇f=(2,2) → next=(0,0)
        r = runner.invoke(main, ['newton', 'x**2 + y**2', '-s', 'x', '1.0', '-s', 'y', '1.0'])
        assert r.exit_code == 0, r.output
        assert '0' in r.output

    def test_pretty(self, runner):
        r = runner.invoke(main, ['newton', '--pretty', 'x**2 + y**2', '-s', 'x', '1.0', '-s', 'y', '1.0'])
        assert r.exit_code == 0, r.output


class TestSucchalv:
    def test_basic(self, runner):
        # f=x²+y², start (2,2): needs at least one halving to satisfy f<8
        r = runner.invoke(main, ['succhalv', 'x**2 + y**2', '2', '2'])
        assert r.exit_code == 0, r.output
        assert 'B*' in r.output

    def test_custom_stepsize(self, runner):
        r = runner.invoke(main, ['succhalv', 'x**2 + y**2', '2', '2', '--stepsize', '2'])
        assert r.exit_code == 0, r.output

    def test_max_steps_limit(self, runner):
        r = runner.invoke(main, ['succhalv', 'x**2 + y**2', '2', '2', '--maxsteps', '1'])
        assert r.exit_code == 0, r.output


class TestAitken:
    def test_basic(self, runner):
        r = runner.invoke(main, ['aitken', '1', '0.5', '0.25'])
        assert r.exit_code == 0, r.output
        assert 'Aitken' in r.output

    def test_longer_sequence(self, runner):
        r = runner.invoke(main, ['aitken', '1', '0.5', '0.25', '0.125'])
        assert r.exit_code == 0, r.output

    def test_too_few_values(self, runner):
        r = runner.invoke(main, ['aitken', '1', '0.5'])
        assert r.exit_code == 0, r.output
        assert 'requires at least 3' in r.output


class TestBroyden:
    def test_basic(self, runner):
        r = runner.invoke(main, ['broyden', 'x**2 + y**2', '1', '1'])
        assert r.exit_code == 0, r.output

    def test_single_step(self, runner):
        r = runner.invoke(main, ['broyden', 'x**2 + y**2', '1', '1', '--steps', '1'])
        assert r.exit_code == 0, r.output

    def test_rational_flag(self, runner):
        r = runner.invoke(main, ['broyden', 'x**2 + y**2', '1', '1', '--rational'])
        assert r.exit_code == 0, r.output

    def test_pretty(self, runner):
        r = runner.invoke(main, ['broyden', '--pretty', 'x**2 + y**2', '1', '1', '--steps', '1'])
        assert r.exit_code == 0, r.output


class TestBroydeninter:
    def test_one_step(self, runner):
        # f=x²+y², start (1,1), ∇f=(2,2), H⁻¹=0.5I
        r = runner.invoke(main, [
            'broydeninter',
            '--startingpoint', '1.0', '1.0',
            '--gradient', '2.0', '2.0',
            '--hessian_inv', '0.5', '0.0', '0.0', '0.5',
            '--steps', '1',
        ])
        assert r.exit_code == 0, r.output

    def test_two_steps(self, runner):
        r = runner.invoke(main, [
            'broydeninter',
            '--startingpoint', '1.0', '1.0',
            '--gradient', '2.0', '2.0',
            '--gradient', '0.0', '0.0',
            '--hessian_inv', '0.5', '0.0', '0.0', '0.5',
            '--steps', '2',
        ])
        assert r.exit_code == 0, r.output

    def test_rational_pretty(self, runner):
        r = runner.invoke(main, [
            'broydeninter',
            '--startingpoint', '1.0', '1.0',
            '--gradient', '2.0', '2.0',
            '--hessian_inv', '0.5', '0.0', '0.0', '0.5',
            '--steps', '1',
            '--rational', '--pretty',
        ])
        assert r.exit_code == 0, r.output


# ── Part 2b ───────────────────────────────────────────────────────────────────

class TestDrawgraph:
    def test_html_output(self, runner):
        with runner.isolated_filesystem():
            r = runner.invoke(main, ['drawgraph', ADJ_CSV])
            assert r.exit_code == 0, r.output
            assert 'graph.html' in r.output

    def test_directed(self, runner):
        with runner.isolated_filesystem():
            r = runner.invoke(main, ['drawgraph', ADJ_CSV, '--directed'])
            assert r.exit_code == 0, r.output

    def test_png_format(self, runner):
        with runner.isolated_filesystem():
            r = runner.invoke(main, ['drawgraph', ADJ_CSV, '--format', 'png'])
            assert r.exit_code == 0, r.output
            assert 'graph.png' in r.output


class TestMst:
    def test_basic(self, runner):
        r = runner.invoke(main, ['mst', ADJ_CSV])
        assert r.exit_code == 0, r.output
        assert 'SUM' in r.output

    def test_total_weight(self, runner):
        # MST of the test graph: C-D:1, B-C:2, A-B:3 → total 6
        r = runner.invoke(main, ['mst', ADJ_CSV])
        assert r.exit_code == 0, r.output
        assert '6' in r.output


class TestDijkstra:
    def test_from_A(self, runner):
        r = runner.invoke(main, ['dijkstra', ADJ_CSV, 'A'])
        assert r.exit_code == 0, r.output
        assert 'Shortest Path' in r.output

    def test_directed(self, runner):
        r = runner.invoke(main, ['dijkstra', ADJ_CSV, 'A', '--directed'])
        assert r.exit_code == 0, r.output

    def test_shortest_path_to_D(self, runner):
        # A→B→C→D costs 3+2+1=6, shorter than direct A→D=7
        r = runner.invoke(main, ['dijkstra', ADJ_CSV, 'A'])
        assert r.exit_code == 0, r.output
        assert '6' in r.output


class TestTraverse:
    def test_bfs(self, runner):
        r = runner.invoke(main, ['traverse', ADJ_CSV, 'A', 'bf'])
        assert r.exit_code == 0, r.output
        assert 'Encounter Order' in r.output

    def test_dfs(self, runner):
        r = runner.invoke(main, ['traverse', ADJ_CSV, 'A', 'df'])
        assert r.exit_code == 0, r.output
        assert 'Encounter Order' in r.output

    def test_bfs_visits_all_nodes(self, runner):
        r = runner.invoke(main, ['traverse', ADJ_CSV, 'A', 'bf'])
        assert r.exit_code == 0, r.output
        for node in ('A', 'B', 'C', 'D'):
            assert node in r.output


class TestFloydwarshall:
    def test_all_pairs(self, runner):
        r = runner.invoke(main, ['floydwarshall', ADJ_CSV])
        assert r.exit_code == 0, r.output

    def test_constrained_nodes(self, runner):
        r = runner.invoke(main, ['floydwarshall', ADJ_CSV, '--onlyuse', 'B, C'])
        assert r.exit_code == 0, r.output


class TestMaxflow:
    def test_basic(self, runner):
        r = runner.invoke(main, ['maxflow', EDGE_CSV, 's', 't'])
        assert r.exit_code == 0, r.output
        assert 'max flow' in r.output

    def test_flow_value(self, runner):
        # s→A:10, s→B:8, A→B:3, A→t:5, B→t:7 → max flow = 12
        r = runner.invoke(main, ['maxflow', EDGE_CSV, 's', 't'])
        assert r.exit_code == 0, r.output
        assert '12' in r.output


class TestMincut:
    def test_edge_list(self, runner):
        r = runner.invoke(main, ['mincut', EDGE_CSV, 's', 't'])
        assert r.exit_code == 0, r.output
        assert 'cut value' in r.output

    def test_adjacency_flag(self, runner):
        r = runner.invoke(main, ['mincut', ADJ_CSV, 'A', 'D', '--adjacency'])
        assert r.exit_code == 0, r.output
        assert 'cut value' in r.output

    def test_cut_value_equals_max_flow(self, runner):
        # by max-flow min-cut theorem, min cut = max flow = 12
        r = runner.invoke(main, ['mincut', EDGE_CSV, 's', 't'])
        assert r.exit_code == 0, r.output
        assert '12' in r.output


class TestMaxmatch:
    def test_basic(self, runner):
        r = runner.invoke(main, ['maxmatch', ADJ_CSV])
        assert r.exit_code == 0, r.output
        assert 'maximum matches' in r.output


class TestMincostmaxflow:
    def test_basic(self, runner):
        r = runner.invoke(main, ['mincostmaxflow', EDGE_CSV, 's', 't'])
        assert r.exit_code == 0, r.output
        assert 'min cost' in r.output
        assert 'max flow' in r.output
