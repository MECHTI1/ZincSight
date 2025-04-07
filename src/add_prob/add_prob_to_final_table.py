from joblib import load
import numpy as np
from src.settings import get_db_connection, KEEP_TEMP_TABLES
from pathlib import Path
import sys

# Environment detection for Colab/Terminal compatibility
try:
    from IPython import get_ipython
    is_notebook = get_ipython() is not None
except:
    is_notebook = False

# Load the Platt scaling model
SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_PATH = SCRIPT_DIR / "platt_model_score_to_prob.joblib"
PLATT_MODEL = load(MODEL_PATH)

def cleanup_tables(conn, cur):
    """Keep only required tables, delete others"""
    required_tables = {
        'minimized_training_cluster_information',
        'motif_representative_coordinates_table_v2',
        'motif_representative_detailed_coordinates_table_v2',
        'training_representative_metal_sites_kruskal_v2'
    }
    cur.execute("SELECT table_name FROM information_schema.tables WHERE table_schema = 'public'")
    for table in cur.fetchall():
        if table[0] not in required_tables:
            cur.execute(f'DROP TABLE IF EXISTS "{table[0]}" CASCADE')
    conn.commit()

def print_bold_message_no_predicted_site_and_cleanup_created_tables():
    """Display bold message and cleanup tables"""
    message = "No predicted zinc-binding sites within the given query structures!"
    print(message, file=sys.stderr, flush=True)  # This ensures visibility
    conn = get_db_connection()
    cur = conn.cursor()
    if not KEEP_TEMP_TABLES:
        cleanup_tables(conn, cur)
    cur.close()
    conn.close()
    return False # Not exist final table -no any predicted sites

def load_scores_to_prob_model_and_predict(scores):
    """Predict probabilities using Platt model"""
    try:
        scores_reshaped = np.array(scores).reshape(-1, 1)
        return PLATT_MODEL.predict_proba(scores_reshaped)[:, 1]
    except Exception as e:
        print(f"Error in probability prediction: {str(e)}")
        return None

def add_column_with_probs():
    """Add probability column to database table"""
    conn = None
    cur = None
    conn = get_db_connection()
    cur = conn.cursor()

    cur.execute("""
        SELECT EXISTS (
            SELECT FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = 'final_compressed_table_with_scored_binding_sites'
        );
        """)
    if not cur.fetchone()[0]:
        print_bold_message_no_predicted_site_and_cleanup_created_tables()

    cur.execute("""
        ALTER TABLE final_compressed_table_with_scored_binding_sites 
        ADD COLUMN IF NOT EXISTS prob REAL;
    """)
    conn.commit()

    cur.execute("SELECT id, score FROM final_compressed_table_with_scored_binding_sites;")
    rows = cur.fetchall()

    if not rows:
        print_bold_message_no_predicted_site_and_cleanup_created_tables()

    ids, scores = zip(*rows)
    probabilities = load_scores_to_prob_model_and_predict(scores)

    if probabilities is None:
        print_bold_message_no_predicted_site_and_cleanup_created_tables()

    updates = [(float(prob), float(id)) for prob, id in zip(probabilities, ids)]
    cur.executemany(
        "UPDATE final_compressed_table_with_scored_binding_sites SET prob = %s WHERE id = %s;",
        updates
    )
    conn.commit()
    return True


if __name__ == "__main__":
    try:
        test_scores = [i for i in range(100)]
        probs = load_scores_to_prob_model_and_predict(test_scores)
        if probs is not None:
            for score, prob in zip(test_scores, probs):
                print(f"Score: {score}, Calibrated Probability: {prob:.4f}")
        add_column_with_probs()
    except Exception as e:
        print_bold_message_no_predicted_site_and_cleanup_created_tables()