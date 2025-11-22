import streamlit as st
import base64
from descriptors_app import show_descriptors_page
from mopac_app import show_mopac_page

# --- Configuración de la página ---
st.set_page_config(
    page_title="ChemY",
    page_icon="assets/logo.png",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- CSS personalizado ---
sidebar_logo_css = """
<style>
.sidebar .sidebar-content {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: flex-start;
}

.sidebar-header {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    width: 100%;
    padding: 1rem 0;
    border-bottom: 1px solid rgba(255, 255, 255, 0.1);
    margin-bottom: 1rem;
}

.sidebar-logo-container {
    width: 90px;
    height: 90px;
    border-radius: 50%;
    background-color: white;
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 5px;
    margin-bottom: 0.5rem;
}

.sidebar-logo {
    width: 100%;
    height: 100%;
    border-radius: 50%;
    object-fit: contain;
}
</style>
"""

st.markdown(sidebar_logo_css, unsafe_allow_html=True)

# --- Mostrar logo con HTML + base64 ---
logo_path = "assets/logo.png"
try:
    with open(logo_path, "rb") as f:
        image_bytes = f.read()
    base64_image = base64.b64encode(image_bytes).decode("utf-8")

    st.sidebar.markdown(
        f"""
        <div class="sidebar-header">
            <div class="sidebar-logo-container">
                <img src="data:image/png;base64,{base64_image}" class="sidebar-logo" alt="ChemY">
            </div>
            <div class="sidebar-title">ChemY</div>
        </div>
        """,
        unsafe_allow_html=True
    )
except FileNotFoundError:
    st.sidebar.warning(f"⚠️ Logo no encontrado en: {logo_path}")
    st.sidebar.markdown('<h2 style="text-align: center; color: white;">ChemY</h2>', unsafe_allow_html=True)

# --- Navegación ---
page = st.sidebar.selectbox("Selecciona una página", ["Visualizador de Descriptores", "MOPAC"])
if page == "Visualizador de Descriptores":
    show_descriptors_page()
elif page == "MOPAC":
    show_mopac_page()