import streamlit as st
import tempfile
import subprocess
import os
import re

# --- Funciones auxiliares para MOPAC (flujo original) ---

def run_mopac(mop_file_path, mopac_executable="MOPAC6J.exe"):
    try:
        st.write(f"🔍 Intentando ejecutar: {mopac_executable} {mop_file_path}")
        if not os.path.exists(mop_file_path):
            st.error(f"❌ El archivo .mop temporal no existe: {mop_file_path}")
            return -1, "", "Archivo .mop no encontrado"

        result = subprocess.run([mopac_executable, mop_file_path], capture_output=True, text=True, timeout=300)
        st.write(f"✅ MOPAC ejecutado. Código de retorno: {result.returncode}")
        return result.returncode, result.stdout, result.stderr
    except FileNotFoundError:
        st.error(f"❌ No se encontró el ejecutable {mopac_executable}. Asegúrate de que esté instalado y accesible.")
        return -1, "", "Executable not found"
    except subprocess.TimeoutExpired:
        st.error("❌ La ejecución de MOPAC tardó demasiado y se interrumpió.")
        return -1, "", "Timeout"
    except Exception as e:
        st.error(f"❌ Error inesperado: {e}")
        return -1, "", str(e)

# --- Nueva función: parsear un .out para comparación ---
def parse_out_for_comparison(out_content, filename):
    data = {
        "name": filename.replace(".out", ""),
        "homo": None,
        "charges": {}
    }

    # Extraer HOMO
    filled_match = re.search(r"NO\. OF FILLED LEVELS\s*=\s*(\d+)", out_content)
    if filled_match:
        num_filled = int(filled_match.group(1))
        eigen_blocks = re.findall(r"ROOT NO\.\s+.*?\n((?:\s*[-+]?\d+\.\d+\s+){1,6})", out_content, re.DOTALL)
        all_energies = []
        for block in eigen_blocks:
            energies = list(map(float, block.split()))
            all_energies.extend(energies)
        if all_energies and (num_filled - 1) < len(all_energies):
            data["homo"] = round(all_energies[num_filled - 1], 5)

    # Extraer cargas netas
    charge_match = re.search(
        r"NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS\s*\n"
        r".*?\n"
        r"((?:\s+\d+\s+[A-Z]\s+[-+]?\d*\.\d+\s+[-+]?\d*\.\d+\s*\n)+)",
        out_content,
        re.DOTALL
    )
    if charge_match:
        charge_lines = charge_match.group(1).strip().split("\n")
        carbon_count = 0
        for line in charge_lines:
            parts = line.split()
            if len(parts) >= 3 and parts[1] == "C":
                carbon_count += 1
                try:
                    charge = float(parts[2])
                    data["charges"][carbon_count] = charge
                except:
                    pass
        if carbon_count != 6:
            data["charges"] = {}

    return data

# --- Nueva función: análisis comparativo de múltiples .out ---
def analyze_multiple_out_files(uploaded_files):
    all_data = []
    benceno_homo = None

    for uploaded_file in uploaded_files:
        content = uploaded_file.getvalue().decode("utf-8")
        parsed = parse_out_for_comparison(content, uploaded_file.name)
        if parsed["homo"] is not None and len(parsed["charges"]) == 6:
            all_data.append(parsed)
            if "benceno" in parsed["name"].lower() and "anilina" not in parsed["name"].lower() and "nitro" not in parsed["name"].lower() and "cl" not in parsed["name"].lower():
                benceno_homo = parsed["homo"]
        else:
            st.warning(f"⚠️ '{uploaded_file.name}' omitido (datos incompletos).")

    if not all_data:
        st.error("❌ No hay datos válidos.")
        return

    # Tabla de HOMO
    st.subheader("📊 Energías del HOMO (eV)")
    homo_table = []
    for d in all_data:
        comparacion = ""
        if benceno_homo is not None and d["homo"] != benceno_homo:
            comparacion = "✅ Activante" if d["homo"] > benceno_homo else "❌ Desactivante"
        homo_table.append({
            "Molécula": d["name"],
            "E(HOMO) [eV]": d["homo"],
            "Clasificación": comparacion
        })
    st.dataframe(homo_table, use_container_width=True)

    # Mostrar cargas POR ÁTOMO (sin asumir posiciones)
    st.subheader("🧬 Cargas Atómicas Netas (C1 a C6)")
    st.markdown("""
    > **Instrucciones para interpretar**:  
    > - Identifica **C1** como el carbono unido al sustituyente.  
    > - Entonces: **C2 y C6 = orto**, **C3 y C5 = meta**, **C4 = para**.  
    > - Compara las cargas en esas posiciones.
    """)
    for d in all_data:
        st.markdown(f"**• {d['name']}**")
        charges = [round(d["charges"][i], 5) for i in range(1, 7)]
        st.text(f"   C1      C2      C3      C4      C5      C6")
        st.text(f"  {charges[0]:7.5f} {charges[1]:7.5f} {charges[2]:7.5f} {charges[3]:7.5f} {charges[4]:7.5f} {charges[5]:7.5f}")
        st.markdown("---")

    if benceno_homo is not None:
        st.info(f"🔹 **Benceno de referencia**: E(HOMO) = {benceno_homo} eV")
# --- Interfaz principal ---
def show_mopac_page():
    st.set_page_config(page_title="ChemY - Análisis de Reactividad Aromática", layout="wide")
    st.title("🧪 ChemY: Análisis Comparativo de Bencenos Monosustituidos")
    st.markdown("**Creado por:** Antonio Elias Sánchez Soto")

    # === SECCIÓN 1: Análisis de múltiples .out ===
    st.header("📤 Comparación de Resultados (.out)")
    st.markdown("""
    Sube los archivos de salida de MOPAC (`.out`) para **benceno, anilina, clorobenceno y nitrobenceno**.
    La app generará tablas comparativas de **E(HOMO)** y **cargas netas** para evaluar reactividad y orientación.
    """)
    uploaded_out_files = st.file_uploader(
        "Sube múltiples archivos .out",
        type=["out"],
        accept_multiple_files=True,
        key="multi_out"
    )
    if uploaded_out_files:
        if len(uploaded_out_files) > 10:
            st.warning("⚠️ Máximo 10 archivos. Se procesarán los primeros 10.")
            uploaded_out_files = uploaded_out_files[:10]
        analyze_multiple_out_files(uploaded_out_files)

    st.markdown("---")

    # === SECCIÓN 2: Generar y ejecutar .mop (flujo original) ===
    st.header("⚙️ Generar y Ejecutar Nuevo Cálculo (.zmt → .mop → .out)")
    with st.expander("📖 Guía de Usuario", expanded=False):
        st.markdown("""
        1. Sube un archivo `.zmt` (geometría de HyperChem).
        2. Selecciona método (PM3, AM1, RM1) y keywords.
        3. Genera el `.mop` y opcionalmente ejecútalo con MOPAC6J.exe.
        """)

    zmt_file = st.file_uploader("Sube archivo .zmt", type=["zmt"])
    if zmt_file is not None:
        st.success("✅ Archivo .zmt cargado.")
        method = st.selectbox("Método", ["PM3", "AM1", "RM1"])
        keywords_input = st.text_input("Keywords adicionales", value="PRECISE VECTORS")
        keywords = f"{method} {keywords_input}".strip()

        if st.button("⚙️ Generar .mop"):
            zmt_content = zmt_file.getvalue().decode("utf-8")
            mop_content = f"{keywords}\n\n" + zmt_content
            st.download_button(
                "📥 Descargar .mop",
                mop_content.encode("utf-8"),
                "entrada.mop",
                "chemical/x-mopac-input"
            )

            with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".mop") as f:
                f.write(mop_content)
                temp_mop = f.name

            if st.button("▶️ Ejecutar MOPAC6J"):
                st.write("⏳ Ejecutando MOPAC...")
                rc, stdout, stderr = run_mopac(temp_mop)
                if rc == 0:
                    out_path = temp_mop.replace(".mop", ".out")
                    if os.path.exists(out_path):
                        with open(out_path, "r") as f:
                            out_content = f.read()
                        st.success("✅ Cálculo completado. Resultados:")
                        st.download_button("📄 Descargar .out", out_content, "resultado.out", "text/plain")

                        # Parseo básico para mostrar
                        parsed = parse_out_for_comparison(out_content, "resultado")
                        if parsed["homo"]:
                            st.write(f"**E(HOMO):** {parsed['homo']} eV")
                        if parsed["charges"]:
                            st.write("**Cargas netas (C1–C6):**", list(parsed["charges"].values()))
                    else:
                        st.error("❌ No se generó el archivo .out.")
                else:
                    st.error(f"❌ Falló MOPAC (código {rc})")
                    if stderr:
                        st.text(stderr)
