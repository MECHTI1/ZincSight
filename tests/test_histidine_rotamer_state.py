"""Rotamer comparisons must not retain exclusions from previous candidates."""
import numpy as np

from src.dis_angles_calculations.with_rotations.return_nonredundant_dicts_his_comb import main


def motif(x):
    return {'A_1': ['HIS', {'ND1': np.array([x, 0., 0.])}]}


def test_previous_candidate_does_not_remove_a_new_conformer():
    assert main([motif(0.)], [motif(0.)]) == []
    distinct = motif(5.)
    result = main([distinct], [motif(0.)])
    assert len(result) == 1
    assert result[0] is distinct


def test_previous_larger_candidate_does_not_leak_indices():
    assert main([motif(0.), motif(0.), motif(0.)], [motif(0.)]) == []
    distinct = motif(5.)
    assert main([distinct], [motif(0.)]) == [distinct]


def test_each_histidine_must_match_for_redundancy():
    candidate = motif(0.)
    candidate['A_2'] = ['HIS', {'ND1': np.array([5., 0., 0.])}]
    previous = motif(0.)
    previous['A_2'] = ['HIS', {'ND1': np.array([0., 0., 0.])}]
    result = main([candidate], [previous])
    assert len(result) == 1
    assert result[0] is candidate
