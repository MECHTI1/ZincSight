"""Regression checks: empty predictions are results; scoring failures are errors."""
from unittest.mock import Mock
import pytest
from src.add_prob import add_prob_to_final_table as scoring


@pytest.mark.parametrize('table_exists,rows', [(False, []), (True, [])])
def test_empty_predictions_return_false(monkeypatch, table_exists, rows):
    conn = Mock()
    cursor = conn.cursor.return_value
    cursor.fetchone.return_value = (table_exists,)
    cursor.fetchall.return_value = rows
    cleanup = Mock(return_value=False)
    monkeypatch.setattr(scoring, 'print_bold_message_no_predicted_site_and_cleanup_created_tables', cleanup)
    assert scoring.add_column_with_probs(conn) is False
    cleanup.assert_called_once()
    cursor.executemany.assert_not_called()
    conn.commit.assert_called()  # Release read locks before cleanup uses another connection.


def test_scoring_failure_is_not_empty_result(monkeypatch):
    conn = Mock()
    cursor = conn.cursor.return_value
    cursor.fetchone.return_value = (True,)
    cursor.fetchall.return_value = [(1, 0.5)]
    monkeypatch.setattr(scoring, 'load_scores_to_prob_model_and_predict', lambda _: None)
    cleanup = Mock(return_value=False)
    monkeypatch.setattr(scoring, 'print_bold_message_no_predicted_site_and_cleanup_created_tables', cleanup)
    with pytest.raises(RuntimeError, match='probability'):
        scoring.add_column_with_probs(conn)
    cleanup.assert_not_called()
