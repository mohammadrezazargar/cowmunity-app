import streamlit as st

from Cowmunity import run_model


TREATMENT_OPTIONS = {
    "None": "no",
    "Imidazole": "imidazole",
    "L-Carnitine": "l-carnitine",
    "Methyl Jasmonate": "methyl jasmonate",
    "Propylpyrazine": "propylpyrazine",
}


st.set_page_config(page_title="Cowmunity Model", layout="centered")
st.title("Cowmunity Model")
st.write("Configure species counts and treatment, then run the model.")

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
