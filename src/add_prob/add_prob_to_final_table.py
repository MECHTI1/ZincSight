from joblib import load
import numpy as np
from src.settings import get_db_connection
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

def print_bold_message():
    """Display bold message in both Colab and terminal"""
    message = "No predicted zinc-binding sites within the given query structures!"
    try:
        if is_notebook:
            from IPython.display import HTML, display
            display(HTML(f"<b>{message}</b>"))
        else:
            print('\033[1m' + message + '\033[0m')
    except:
        print(message)  # Fallback
    finally:
        sys.exit(0)

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
    try:
        conn = get_db_connection()
        cur = conn.cursor()

        # Check if table exists
        cur.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = 'public' 
                AND table_name = 'final_compressed_table_with_scored_binding_sites'
            );
        """)
        if not cur.fetchone()[0]:
            print_bold_message()

        # Add probability column
        cur.execute("""
            ALTER TABLE final_compressed_table_with_scored_binding_sites 
            ADD COLUMN IF NOT EXISTS prob REAL;
        """)
        conn.commit()

        # Get existing data
        cur.execute("SELECT id, score FROM final_compressed_table_with_scored_binding_sites;")
        rows = cur.fetchall()

        if not rows:  # Empty table check
            print_bold_message()

        # Calculate probabilities
        ids, scores = zip(*rows)
        probabilities = load_scores_to_prob_model_and_predict(scores)

        if probabilities is None:
            print_bold_message()

        # Update probabilities
        updates = [(prob, id) for prob, id in zip(probabilities, ids)]
        cur.executemany(
            "UPDATE final_compressed_table_with_scored_binding_sites SET prob = %s WHERE id = %s;",
            updates
        )
        conn.commit()

    except Exception as e:
        print_bold_message()
    finally:
        if cur:
            try:
                cur.close()
            except:
                pass
        if conn:
            try:
                conn.close()
            except:
                pass

if __name__ == "__main__":
    try:
        # Test section
        test_scores = [i for i in range(100)]
        probs = load_scores_to_prob_model_and_predict(test_scores)
        if probs is not None:
            for score, prob in zip(test_scores, probs):
                print(f"Score: {score}, Calibrated Probability: {prob:.4f}")

        add_column_with_probs()
    except Exception as e:
        print_bold_message()