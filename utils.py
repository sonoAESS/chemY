from rdkit import Chem

def load_molecule_from_file(file_buffer):
    """
    Carga una molécula desde un archivo subido en Streamlit.
    """
    try:
        with open("temp_mol.mol", "w") as f:
            f.write(file_buffer.getvalue().decode("utf-8"))
        mol = Chem.MolFromMolFile("temp_mol.mol", removeHs=False, sanitize=True)
        return mol
    except Exception:
        return None