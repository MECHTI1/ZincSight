from joblib import load
import numpy as np
from src.settings import get_db_connection
from pathlib import Path

# Load the Platt scaling model once at the start for efficiency
SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_PATH = SCRIPT_DIR / "platt_model_score_to_prob.joblib"
PLATT_MODEL = load(MODEL_PATH)

def load_scores_to_prob_model_and_predict(scores):
    """
    Predict calibrated probabilities for a list of scores using the loaded Platt model.
    """
    # Ensure scores is a 2D array (n, 1) for the model
    scores_reshaped = np.array(scores).reshape(-1, 1)

    # Predict calibrated probabilities
    calibrated_probs = PLATT_MODEL.predict_proba(scores_reshaped)[:, 1]

    return calibrated_probs


def add_column_with_probs():
    """
    Add a new column to the database table and populate it with probabilities.
    """
    # Get database connection
    conn = get_db_connection()
    cur = conn.cursor()

    # Add a new column for probabilities if it doesn't exist
    cur.execute("""
        ALTER TABLE final_compressed_table_with_scored_binding_sites 
        ADD COLUMN IF NOT EXISTS prob REAL;
    """)
    conn.commit()

    # Select scores and IDs from the table
    cur.execute("SELECT id, score FROM final_compressed_table_with_scored_binding_sites;")
    rows = cur.fetchall()

    # Process scores in bulk
    ids, scores = zip(*rows)  # Unzip rows into separate lists
    probabilities = load_scores_to_prob_model_and_predict(scores)

    # Prepare update statements
    updates = [(prob, id) for prob, id in zip(probabilities, ids)]

    # Update the table in bulk
    cur.executemany(
        "UPDATE final_compressed_table_with_scored_binding_sites SET prob = %s WHERE id = %s;",
        updates
    )
    conn.commit()

if __name__ == "__main__":

    # Test multiple scores
    test_scores = [i for i in range(100)]
    probs = load_scores_to_prob_model_and_predict(test_scores)
    for score, prob in zip(test_scores, probs):
        print(f"Score: {score}, Calibrated Probability: {prob:.4f}")

    # Add probabilities to the database
    add_column_with_probs()
