import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
import tempfile
import re
from collections import defaultdict

# Importar RDKit con manejo de errores
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError as e:
    st.error(f"Error importing RDKit: {e}")
    RDKIT_AVAILABLE = False

# Importar Draw por separado con manejo de errores
try:
    from rdkit.Chem import Draw
    DRAW_AVAILABLE = True
except ImportError:
    DRAW_AVAILABLE = False

# Importar py3Dmol con manejo de errores
try:
    import py3Dmol
    PY3DMOL_AVAILABLE = True
except ImportError:
    PY3DMOL_AVAILABLE = False

def mol_to_3d_view(mol, width=600, height=600):
    """Convierte molécula RDKit a visualización 3D con py3Dmol"""
    if not PY3DMOL_AVAILABLE or mol is None:
        return None
    
    try:
        mol_with_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDGv3())
        mol_block = Chem.MolToMolBlock(mol_with_h)
        view = py3Dmol.view(width=width, height=height)
        view.addModel(mol_block, "mol")
        view.setStyle({"stick": {}})
        view.setBackgroundColor("white")
        view.zoomTo()
        return view
    except Exception as e:
        st.error(f"Error generating 3D view: {e}")
        return None

# 🔬 NUEVA FUNCIÓN: Parsear archivo .out de MOPAC (PM3, AM1, RM1 con VECTORS)
def parse_mopac_pm3_out(content_lines):
    charges = {}
    in_charges = False
    for line in content_lines:
        if "NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS" in line:
            in_charges = True
            continue
        if in_charges and "ATOM NO." in line:
            continue
        if in_charges and line.strip() == "":
            break
        if in_charges:
            match = re.match(r'\s*(\d+)\s+([A-Z])\s+([+-]?\d*\.\d+)', line)
            if match:
                atom_idx = int(match.group(1))
                elem = match.group(2)
                charge = float(match.group(3))
                if elem == "C":
                    charges[atom_idx] = charge

    # Extraer eigenvalues
    eigenvalues = []
    in_eigen = False
    for line in content_lines:
        if "ROOT NO." in line:
            in_eigen = True
            continue
        if in_eigen and line.strip() == "":
            in_eigen = False
            continue
        if in_eigen:
            # Solo líneas con números en notación decimal
            parts = line.split()
            if parts and all(re.match(r'[+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?', p) for p in parts if p):
                try:
                    vals = [float(p) for p in parts]
                    eigenvalues.extend(vals)
                except ValueError:
                    continue

    # Determinar HOMO y pre-HOMO (asumimos sistema cerrado, 12 orbitales ocupados para benceno)
    homo_energy = prehomo_energy = None
    if len(eigenvalues) >= 12:
        homo_energy = eigenvalues[11]      # ROOT 12
        prehomo_energy = eigenvalues[10]   # ROOT 11

    # Extraer coeficientes pz del bloque con ROOT 7-12
    pz_coeffs = {}
    i = 0
    while i < len(content_lines):
        line = content_lines[i]
        if "ROOT NO." in line:
            root_nums = list(map(int, re.findall(r'\d+', line)))
            if 7 in root_nums and 12 in root_nums:
                # Saltar la línea de energía
                i += 2
                atom_data = {}
                while i < len(content_lines) and content_lines[i].strip() != "":
                    l = content_lines[i]
                    if "PZ" in l:
                        parts = l.split()
                        if len(parts) >= 9:  # "PZ C idx + 6 coeficientes"
                            try:
                                atom_idx = int(parts[2])
                                # Columnas: ROOT7→idx3, ROOT8→4, ..., ROOT12→8
                                prehomo_coef = float(parts[7])  # ROOT 11 → columna 5 del bloque (índice 7)
                                homo_coef = float(parts[8])     # ROOT 12 → columna 6 (índice 8)
                                atom_data[atom_idx] = {
                                    "prehomo_coef": prehomo_coef,
                                    "homo_coef": homo_coef,
                                    "prehomo_sq": prehomo_coef**2,
                                    "homo_sq": homo_coef**2
                                }
                            except (ValueError, IndexError):
                                pass
                    i += 1
                pz_coeffs = atom_data
                break
            else:
                # Saltar bloque
                i += 2
                while i < len(content_lines) and content_lines[i].strip() != "":
                    i += 1
        i += 1

    return {
        "charges": charges,
        "homo_energy": homo_energy,
        "prehomo_energy": prehomo_energy,
        "pz_coeffs": pz_coeffs
    }

def show_descriptors_page():
    st.title("🔬 Visualizador de Descriptores Moleculares")
    st.markdown("**Creado por:** Antonio Elias Sánchez Soto")

    if not RDKIT_AVAILABLE:
        st.error("❌ RDKit no está disponible en este entorno.")
        return

    # --- Guía de usuario (sin cambios) ---
    with st.expander("📖 Guía de Usuario", expanded=False):
        st.markdown("""
        ### 📥 Cómo cargar archivos
        ... (tu guía original permanece igual) ...
        """)

    # --- Sidebar ---
    st.sidebar.header("📁 Subir archivos")

    # Archivo .mol
    mol_file = st.sidebar.file_uploader("Sube tu archivo .mol", type=["mol"])
    mol = None
    if mol_file is not None:
        try:
            mol_data = mol_file.getvalue().decode("utf-8")
            mol = Chem.MolFromMolBlock(mol_data, removeHs=False, sanitize=True)
            if mol is None:
                st.sidebar.error("❌ No se pudo cargar la molécula.")
            else:
                st.sidebar.success("✅ Molécula cargada correctamente.")
        except Exception as e:
            st.sidebar.error(f"❌ Error al cargar la molécula: {e}")
            mol = None
    else:
        st.sidebar.warning("⚠️ Sube un archivo .mol")

    # Archivos de descriptores
    desc_files = st.sidebar.file_uploader(
        "Sube tus archivos de descriptores (.csv o .txt)", 
        type=["csv", "txt"], 
        accept_multiple_files=True
    )

    # Archivo .out de MOPAC
    mopac_file = st.sidebar.file_uploader("Sube tu archivo de salida de MOPAC (.out)", type=["out"])

    # --- Procesar descriptores (sin cambios) ---
    descriptors_data = {}
    if desc_files:
        for file in desc_files:
            try:
                if file.name.endswith(".txt"):
                    df_temp = pd.read_csv(file, sep="\t", skiprows=2)
                    df_temp.columns = df_temp.columns.str.strip()
                    df_temp.replace(-999, pd.NA, inplace=True)
                    df_temp.dropna(axis=1, how='all', inplace=True)
                    fuente_temp = "e-Dragon"
                elif file.name.endswith(".csv"):
                    df_temp = pd.read_csv(file, index_col=0)
                    df_temp.replace(-999, pd.NA, inplace=True)
                    if 'Name' in df_temp.columns:
                        fuente_temp = "padelpy"
                    else:
                        fuente_temp = "PaDEL"
                else:
                    continue
                descriptors_data[file.name] = {"df": df_temp, "fuente": fuente_temp}
                st.sidebar.success(f"✅ {file.name} ({fuente_temp}) cargado.")
            except Exception as e:
                st.sidebar.error(f"❌ Error al cargar {file.name}: {e}")
                continue

    df = pd.DataFrame()
    fuente = None
    if descriptors_data:
        selected_file = st.sidebar.selectbox(
            "Selecciona archivo de descriptores para visualizar",
            list(descriptors_data.keys())
        )
        df = descriptors_data[selected_file]["df"]
        fuente = descriptors_data[selected_file]["fuente"]

    # --- Procesar archivo MOPAC ---
    mopac_results = None
    if mopac_file is not None:
        try:
            content = mopac_file.getvalue().decode("utf-8")
            lines = content.splitlines()
            mopac_results = parse_mopac_pm3_out(lines)
            st.sidebar.success("✅ Archivo MOPAC (.out) procesado.")
        except Exception as e:
            st.sidebar.error(f"❌ Error al procesar archivo MOPAC: {e}")

    # --- Validar carga ---
    if df.empty and mol is None and mopac_results is None:
        st.warning("⚠️ Sube al menos un archivo (.mol, descriptores o .out de MOPAC) para continuar.")
        return

    # --- Tabs ---
    tab_names = ["📊 Descriptores"]
    if mol is not None:
        tab_names.extend(["🖼️ 2D", "🧬 3D"])
    if mopac_results is not None:
        tab_names.append("🔬 MOPAC Analysis")

    tabs = st.tabs(tab_names)

    # Pestaña Descriptores
    if "📊 Descriptores" in tab_names:
        with tabs[0]:
            if not df.empty:
                mol_id = st.sidebar.selectbox("Selecciona molécula", df.index, key="desc_mol")
                df_filtered = df.loc[[mol_id]].copy()
                desc_search = st.sidebar.text_input("Buscar descriptor (ej. MW, nHDon)", key="desc_search")
                if desc_search:
                    matching_cols = [c for c in df_filtered.columns if desc_search.upper() in c.upper()]
                    df_filtered = df_filtered[matching_cols]

                st.subheader(f"Descriptores ({fuente}) para {mol_id}")
                descriptors_series = df_filtered.T
                descriptors_series.columns = ['Valor']
                st.dataframe(descriptors_series)
                st.download_button(
                    label="📥 Descargar descriptores",
                    data=descriptors_series.to_csv().encode("utf-8"),
                    file_name=f"descriptores_{fuente}_{mol_id}.csv",
                    mime="text/csv",
                )
            else:
                st.info("Sube un archivo de descriptores para ver esta pestaña.")

    # Pestaña 2D
    if "🖼️ 2D" in tab_names:
        idx = tab_names.index("🖼️ 2D")
        with tabs[idx]:
            st.subheader("Estructura 2D de la molécula")
            if mol and DRAW_AVAILABLE:
                try:
                    img = Draw.MolToImage(mol, size=(600, 600))
                    st.image(img, caption="Estructura 2D", use_container_width=True)
                except Exception as e:
                    st.error(f"❌ Error al generar imagen 2D: {e}")
            else:
                st.warning("⚠️ No hay molécula cargada o falta RDKit Draw.")

    # Pestaña 3D
    if "🧬 3D" in tab_names:
        idx = tab_names.index("🧬 3D")
        with tabs[idx]:
            st.subheader("Estructura 3D de la molécula")
            if mol and PY3DMOL_AVAILABLE:
                try:
                    view = mol_to_3d_view(mol, width=600, height=600)
                    if view:
                        html_content = view._make_html()
                        components.html(html_content, width=600, height=600, scrolling=False)
                    else:
                        st.warning("❌ No se pudo generar la visualización 3D.")
                except Exception as e:
                    st.error(f"❌ Error al generar visualización 3D: {e}")
            else:
                st.warning("⚠️ No hay molécula cargada o falta py3Dmol.")

    # 🔬 NUEVA Pestaña: MOPAC Analysis
    if "🔬 MOPAC Analysis" in tab_names:
        idx = tab_names.index("🔬 MOPAC Analysis")
        with tabs[idx]:
            st.subheader("Análisis de Resultados MOPAC")
            if mopac_results:
                # Mostrar energías
                st.markdown(f"**Energía HOMO:** `{mopac_results['homo_energy']:.5f}` eV")
                st.markdown(f"**Energía pre-HOMO:** `{mopac_results['prehomo_energy']:.5f}` eV")
                st.markdown("---")

                # Generar tabla
                rows = []
                for atom_idx in sorted(mopac_results["charges"].keys()):
                    charge = mopac_results["charges"][atom_idx]
                    coef_data = mopac_results["pz_coeffs"].get(atom_idx, {
                        "homo_coef": 0.0, "prehomo_coef": 0.0,
                        "homo_sq": 0.0, "prehomo_sq": 0.0
                    })
                    rows.append({
                        "Átomo (C#)": atom_idx,
                        "Carga neta": f"{charge:.4f}",
                        "Coef. pz (HOMO)": f"{coef_data['homo_coef']:.5f}",
                        "(Coef. pz)² (HOMO)": f"{coef_data['homo_sq']:.5f}",
                        "Coef. pz (pre-HOMO)": f"{coef_data['prehomo_coef']:.5f}",
                        "(Coef. pz)² (pre-HOMO)": f"{coef_data['prehomo_sq']:.5f}"
                    })

                df_mopac = pd.DataFrame(rows)
                st.dataframe(df_mopac, use_container_width=True)

                # Botón de descarga
                st.download_button(
                    label="📥 Descargar análisis MOPAC",
                    data=df_mopac.to_csv(index=False).encode("utf-8"),
                    file_name="mopac_analysis.csv",
                    mime="text/csv"
                )
            else:
                st.info("No se pudo cargar el análisis de MOPAC.")
