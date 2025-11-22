import subprocess
import sys
import os

# Ruta al archivo principal de Streamlit
main_script = "main.py"

if __name__ == "__main__":
    # Llama a Streamlit como si ejecutáramos "streamlit run main.py"
    subprocess.run([sys.executable, "-m", "streamlit", "run", main_script] + sys.argv[1:])