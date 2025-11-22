import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
import tempfile

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
        # Añadir hidrógenos y generar conformación 3D
        mol_with_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDGv3())
        
        # Convertir a formato MOL
        mol_block = Chem.MolToMolBlock(mol_with_h)
        
        # Crear visualización 3D
        view = py3Dmol.view(width=width, height=height)
        view.addModel(mol_block, "mol")
        view.setStyle({"stick": {}})
        view.setBackgroundColor("white")
        view.zoomTo()
        
        return view
    except Exception as e:
        st.error(f"Error generating 3D view: {e}")
        return None

def show_descriptors_page():
    # --- Título ---
    st.title("🔬 Visualizador de Descriptores Moleculares")
    st.markdown("**Creado por:** Antonio Elias Sánchez Soto")

    # Verificar disponibilidad de librerías
    if not RDKIT_AVAILABLE:
        st.error("""
        ❌ RDKit no está disponible en este entorno. 
        Esta aplicación requiere RDKit para funcionar correctamente.
        """)
        return

    # --- Guía de usuario ---
    with st.expander("📖 Guía de Usuario", expanded=False):
        st.markdown("""
        ### 📥 Cómo cargar archivos

        1. **Cargar archivo .mol**:
           - En el panel lateral izquierdo (sidebar), use el selector "Sube tu archivo .mol".
           - Seleccione un archivo con extensión `.mol` que contenga la estructura de la molécula.
           - Este archivo se usará para generar las visualizaciones 2D y 3D.

        2. **Cargar archivos de descriptores**:
           - Use el selector "Sube tus archivos de descriptores (.csv o .txt)".
           - Puede subir uno o más archivos a la vez.
           - Acepta archivos `.csv` (por ejemplo, de PaDEL o padelpy) o `.txt` (por ejemplo, de e-Dragon).
           - La aplicación detectará automáticamente el tipo de software (e-Dragon, PaDEL, padelpy) basándose en el formato del archivo.

        --- 

        ### 🧪 Cómo ver las visualizaciones

        1. **Visualización 2D**:
           - Vaya a la pestaña **🖼️ 2D**.
           - Se mostrará una imagen 2D de la molécula cargada.

        2. **Visualización 3D**:
           - Vaya a la pestaña **🧬 3D**.
           - Se mostrará una vista interactiva 3D de la molécula. Puedes hacer zoom, rotarla, etc.

        ---

        ### 📊 Cómo ver y filtrar los descriptores

        1. **Seleccionar archivo de descriptores**:
           - En el sidebar, después de cargar archivos, aparecerá un selector: "Selecciona archivo de descriptores para visualizar".
           - Elija el archivo del cual desea ver los descriptores.

        2. **Seleccionar molécula**:
           - Si el archivo contiene múltiples moléculas, use el selector "Selecciona molécula".

        3. **Filtrar descriptores**:
           - Use la caja de texto "Buscar descriptor (ej. MW, nHDon)" para filtrar los descriptores mostrados.

        4. **Descargar descriptores**:
           - En la pestaña **📊 Descriptores**, hay un botón "📥 Descargar descriptores" para guardar los datos visualizados en formato CSV.

        """)

    # --- Sidebar ---
    st.sidebar.header("📁 Subir archivos")

    # Subir molécula
    mol_file = st.sidebar.file_uploader("Sube tu archivo .mol", type=["mol"])
    mol = None
    
    if mol_file is not None:
        try:
            # Leer directamente desde el archivo subido
            mol_data = mol_file.getvalue().decode("utf-8")
            mol = Chem.MolFromMolBlock(mol_data, removeHs=False, sanitize=True)
            
            if mol is None:
                # Intentar alternativas si falla
                mol = Chem.MolFromMolFile(mol_file, removeHs=False, sanitize=True)
                
            if mol is None:
                st.sidebar.error("❌ No se pudo cargar la molécula. Verifique el formato del archivo.")
            else:
                st.sidebar.success("✅ Molécula cargada correctamente.")
        except Exception as e:
            st.sidebar.error(f"❌ Error al cargar la molécula: {e}")
            mol = None
    else:
        st.sidebar.warning("⚠️ Sube un archivo .mol")

    # Subir múltiples archivos de descriptores
    desc_files = st.sidebar.file_uploader(
        "Sube tus archivos de descriptores (.csv o .txt)", 
        type=["csv", "txt"], 
        accept_multiple_files=True
    )

    # Diccionario para almacenar DataFrames
    descriptors_data = {}

    if desc_files:
        for file in desc_files:
            try:
                # Detectar tipo de archivo
                if file.name.endswith(".txt"):
                    # Asumir e-Dragon
                    df_temp = pd.read_csv(file, sep="\t", skiprows=2)
                    df_temp.columns = df_temp.columns.str.strip()
                    df_temp.replace(-999, pd.NA, inplace=True)
                    df_temp.dropna(axis=1, how='all', inplace=True)
                    fuente_temp = "e-Dragon"
                elif file.name.endswith(".csv"):
                    # Asumir PaDEL o padelpy
                    df_temp = pd.read_csv(file, index_col=0)
                    df_temp.replace(-999, pd.NA, inplace=True)
                    if 'Name' in df_temp.columns:
                        fuente_temp = "padelpy"
                    else:
                        fuente_temp = "PaDEL"
                else:
                    continue

                # Almacenar en el diccionario
                descriptors_data[file.name] = {"df": df_temp, "fuente": fuente_temp}
                st.sidebar.success(f"✅ {file.name} ({fuente_temp}) cargado.")
                
            except Exception as e:
                st.sidebar.error(f"❌ Error al cargar {file.name}: {e}")
                continue

    # Verificar si se cargó algún archivo de descriptores
    if descriptors_data:
        # Selector de archivo de descriptores
        selected_file = st.sidebar.selectbox(
            "Selecciona archivo de descriptores para visualizar",
            list(descriptors_data.keys())
        )

        # Obtener el DataFrame y la fuente del archivo seleccionado
        df = descriptors_data[selected_file]["df"]
        fuente = descriptors_data[selected_file]["fuente"]
    else:
        df = pd.DataFrame()
        fuente = None
        if desc_files:  # Solo mostrar advertencia si se subieron archivos pero no se pudieron cargar
            st.sidebar.error("❌ No se pudo cargar ningún archivo de descriptores.")

    # --- Validar carga ---
    if df.empty or mol is None:
        st.warning("⚠️ Sube tanto el archivo .mol como al menos un archivo de descriptores para continuar.")
        return

    # --- Sidebar ---
    st.sidebar.header("🔍 Filtros")

    # Seleccionar molécula
    mol_id = st.sidebar.selectbox("Selecciona molécula", df.index)
    df_filtered = df.loc[[mol_id]].copy()

    # Filtrar descriptores
    desc_search = st.sidebar.text_input("Buscar descriptor (ej. MW, nHDon)")
    if desc_search:
        matching_cols = [c for c in df_filtered.columns if desc_search.upper() in c.upper()]
        df_filtered = df_filtered[matching_cols]

    # --- Tabs ---
    tab1, tab2, tab3 = st.tabs(["📊 Descriptores", "🖼️ 2D", "🧬 3D"])

    with tab1:
        st.subheader(f"Descriptores ({fuente}) para {mol_id} (de {selected_file})")
        # Mostrar la fila de descriptores como una tabla vertical (descriptor | valor)
        descriptors_series = df_filtered.T
        descriptors_series.columns = ['Valor']
        st.dataframe(descriptors_series)
        
        # Botón de descarga
        st.download_button(
            label="📥 Descargar descriptores",
            data=descriptors_series.to_csv().encode("utf-8"),
            file_name=f"descriptores_{fuente}_{mol_id}.csv",
            mime="text/csv",
        )

    with tab2:
        st.subheader("Estructura 2D de la molécula")
        if mol and DRAW_AVAILABLE:
            try:
                img = Draw.MolToImage(mol, size=(600, 600))
                st.image(img, caption="Estructura 2D", use_container_width=True)
            except Exception as e:
                st.error(f"❌ Error al generar imagen 2D: {e}")
                st.info("La visualización 2D no está disponible en este entorno.")
        else:
            if not DRAW_AVAILABLE:
                st.warning("🔧 La visualización 2D no está disponible en este entorno.")
            else:
                st.warning("⚠️ No hay molécula cargada para visualizar.")

    with tab3:
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
            if not PY3DMOL_AVAILABLE:
                st.warning("🔧 La visualización 3D no está disponible. Instala py3dmol.")
            else:
                st.warning("⚠️ No hay molécula cargada para visualizar.")

# Información del estado de las librerías
def show_library_status():
    st.sidebar.markdown("---")
    st.sidebar.subheader("🔧 Estado de Librerías")
    
    status_icon = {
        True: "✅",
        False: "❌"
    }
    
    st.sidebar.write(f"{status_icon[RDKIT_AVAILABLE]} RDKit")
    st.sidebar.write(f"{status_icon[DRAW_AVAILABLE]} RDKit Draw")
    st.sidebar.write(f"{status_icon[PY3DMOL_AVAILABLE]} py3Dmol")

# Mostrar estado de librerías
show_library_status()