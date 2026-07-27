import subprocess
import sys

import streamlit as st

st.set_page_config(page_title="Cowmunity Model", layout="centered")


# --- One-time GAMSPy license setup (reads from Streamlit secrets) -----------
# On Streamlit Community Cloud there's no interactive terminal, so the GAMSPy
# license has to be installed programmatically on first run of each new
# container. Add your license under Settings -> Secrets as:
#
#   GAMSPY_LICENSE = "your-36-character-access-code-or-6-line-license-text"
#
@st.cache_resource
def _install_gamspy_license():
    license_value = st.secrets.get("GAMSPY_LICENSE", None)
    if not license_value:
        return "missing", "", ""
    result = subprocess.run(
        [sys.executable, "-m", "gamspy", "install", "license", license_value],
        capture_output=True,
        text=True,
    )
    status = "installed" if result.returncode == 0 else f"error (code {result.returncode})"
    return status, result.stdout, result.stderr


license_status, license_stdout, license_stderr = _install_gamspy_license()

# Always print to the server log, regardless of outcome, so it shows up in
# "Manage app" logs on Streamlit Cloud.
print(f"GAMSPY LICENSE INSTALL STATUS: {license_status}")
print(f"GAMSPY LICENSE INSTALL STDOUT: {license_stdout}")
print(f"GAMSPY LICENSE INSTALL STDERR: {license_stderr}")

# Also show a quick summary of what gamspy currently thinks is installed.
try:
    check = subprocess.run(
        [sys.executable, "-m", "gamspy", "show", "license"],
        capture_output=True,
        text=True,
    )
    print(f"GAMSPY SHOW LICENSE STDOUT: {check.stdout}")
    print(f"GAMSPY SHOW LICENSE STDERR: {check.stderr}")
except Exception as e:
    print(f"Could not run 'gamspy show license': {e}")

if license_status == "missing":
    st.warning(
        "No GAMSPY_LICENSE found in Streamlit secrets. The app will run on the "
        "free demo license, which caps NLP models at 2,500 variables/constraints — "
        "too small for this model with more than one organism copy. "
        "Add your free academic GAMSPy license in the app's Secrets settings."
    )
elif license_status.startswith("error"):
    st.error(f"Could not install GAMSPy license: {license_status}. Check Manage app logs for details.")
# ------------------------------------------------------------------------

from Cowmunity import run_model  # noqa: E402  (import after license setup)


TREATMENT_OPTIONS = {
    "None": "no",
    "Imidazole": "imidazole",
    "L-Carnitine": "l-carnitine",
    "Methyl Jasmonate": "methyl jasmonate",
    "Propylpyrazine": "propylpyrazine",
}


st.title("Cowmunity Model")
st.write("Configure species counts and treatment, then run the model.")

with st.expander("Diagnostics: test GAMSPy with a minimal unrelated model"):
    if st.button("Run minimal GAMSPy test"):
        try:
            import gamspy as gp

            m = gp.Container()
            i = gp.Set(m, "i", records=["seattle", "san-diego"])
            j = gp.Set(m, "j", records=["new-york", "chicago", "topeka"])
            a = gp.Parameter(m, "a", domain=i, records=[("seattle", 350), ("san-diego", 600)])
            b = gp.Parameter(m, "b", domain=j, records=[("new-york", 325), ("chicago", 300), ("topeka", 275)])
            d = gp.Parameter(
                m, "d", domain=[i, j],
                records=[
                    ("seattle", "new-york", 2.5), ("seattle", "chicago", 1.7), ("seattle", "topeka", 1.8),
                    ("san-diego", "new-york", 2.5), ("san-diego", "chicago", 1.8), ("san-diego", "topeka", 1.4),
                ],
            )
            f = gp.Parameter(m, "f", records=90)
            c = gp.Parameter(m, "c", domain=[i, j])
            c[i, j] = f * d[i, j] / 1000
            x = gp.Variable(m, "x", domain=[i, j], type="Positive")
            supply = gp.Equation(m, "supply", domain=i)
            supply[i] = gp.Sum(j, x[i, j]) <= a[i]
            demand = gp.Equation(m, "demand", domain=j)
            demand[j] = gp.Sum(i, x[i, j]) >= b[j]
            cost = gp.Sum((i, j), c[i, j] * x[i, j])
            transport = gp.Model(
                m, name="transport", equations=m.getEquations(),
                problem="LP", sense=gp.Sense.MIN, objective=cost,
            )
            transport.solve()
            st.success(f"Minimal GAMSPy test SOLVED. Objective value: {cost.toValue()}")
        except Exception as e:
            st.error(f"Minimal GAMSPy test FAILED: {e}")
            print(f"MINIMAL GAMSPY TEST FAILED: {e}")

col1, col2, col3 = st.columns(3)
with col1:
    mgk_count = st.number_input("MGK count", min_value=1, value=1, step=1)
with col2:
    prm_count = st.number_input("PRM count", min_value=1, value=1, step=1)
with col3:
    rfl_count = st.number_input("RFL count", min_value=1, value=1, step=1)

treatment_label = st.selectbox("Treatment", list(TREATMENT_OPTIONS.keys()))
variable_choice = st.selectbox("Objective", ["biomass_outer", "ATP_outer"])

run_clicked = st.button("Run model")

if run_clicked:
    species_counts = {"mgk": mgk_count, "prm": prm_count, "rfl": rfl_count}
    treatment = TREATMENT_OPTIONS[treatment_label]
    with st.spinner("Running model..."):
        results, output_dir = run_model(
            species_counts=species_counts,
            treatment=treatment,
            variable_choice=variable_choice,
        )
    st.success("Run completed.")
    st.subheader("Results")
    st.json(results)
    st.subheader("Output files")
    st.code(output_dir)
    st.write("Files saved:", ["mgk_records.csv", "prm_records.csv", "rfl_records.csv"])
