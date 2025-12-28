import streamlit as st
import os
import tempfile
import csv
from Bio import SeqIO
import json
import shutil

from LRS import (
    extract_sequences,
    build_command,
    simulate_badreads,
    combine_and_sample_reads,
    Nanoplot,
    plot_sample_composition_plotly,
)

st.set_page_config(page_title="Long Read Sampler")

st.title("Long Read Sampler (LRS) GUI")

# upload FASTA file
fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa"])

# load default quality presets JSON
default_quality_file = "quality_presets.json"
with open(default_quality_file) as f:
    quality_presets = json.load(f)

# upload custom quality presets JSON
quality_file = st.file_uploader("Upload Quality Presets JSON", type=["json"])
if quality_file:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".json") as tmp_json:
        tmp_json.write(quality_file.read())
        quality_path = tmp_json.name
    with open(quality_path) as f:
        quality_presets = json.load(f)

if fasta_file:
    # get headers from fasta file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp_fasta:
        tmp_fasta.write(fasta_file.read())
        fasta_path = tmp_fasta.name

    all_headers = [seq.description for seq in SeqIO.parse(fasta_path, "fasta")]

    # optional: load a CSV/TSV file that contains species and percentage columns (collapsible)
    with st.expander("Import composition (CSV/TSV)", expanded=False):
        composition_file = st.file_uploader("Upload composition CSV/TSV (species + percentage)", type=["csv", "tsv", "txt"])    
        comp_preview = None
        if composition_file:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp_comp:
                tmp_comp.write(composition_file.read())
                comp_path = tmp_comp.name

            # attempt to detect delimiter and read header
            try:
                with open(comp_path, newline='', encoding='utf-8') as f:
                    sample = f.read(2048)
                    f.seek(0)
                    try:
                        dialect = csv.Sniffer().sniff(sample, delimiters=[',', '\t', ';'])
                    except Exception:
                        dialect = csv.get_dialect('excel')
                        dialect.delimiter = ','
                    reader = csv.DictReader(f, dialect=dialect)
                    fieldnames = reader.fieldnames or []

                if not fieldnames:
                    st.warning('Uploaded file does not contain a header row; column selection requires headers.')
                else:
                    st.markdown('**Import composition from file**')
                    species_col = st.selectbox('Species column', fieldnames, index=0)
                    perc_col = st.selectbox('Percentage column', fieldnames, index=min(1, len(fieldnames)-1))
                    total_reads = st.number_input('Total reads to distribute', min_value=1, value=1000)
                    default_preset_for_file = st.selectbox('Default quality preset for imported species', list(quality_presets.keys()))

                    if st.button('Import composition file'):
                        added = 0
                        with open(comp_path, newline='', encoding='utf-8') as f:
                            reader = csv.DictReader(f, dialect=dialect)
                            for row in reader:
                                name = row.get(species_col, '').strip()
                                perc_raw = row.get(perc_col, '').strip()
                                if not name or not perc_raw:
                                    continue
                                # allow percentages like '12.5' or '12.5%'
                                try:
                                    perc = float(perc_raw.strip().rstrip('%'))
                                except Exception:
                                    continue
                                number = max(1, int(round(total_reads * perc / 100.0)))
                                st.session_state.targets.append({
                                    "header": name,
                                    "number": int(number),
                                    "sim_quality": default_preset_for_file,
                                })
                                added += 1

                        st.success(f'Imported {added} species from composition file.')
            except:
                st.error('Failed to read the uploaded composition file. Please ensure it is a valid CSV/TSV file.')

    st.markdown('#') 

    st.subheader("Select sequences and simulation parameters")

    # session_state to store selected targets
    if "targets" not in st.session_state:
        st.session_state.targets = []

    # select sequences, qualities and quantities for simulation and sampling
    with st.form("add_sequence_form"):
        seq_choice = st.selectbox("Select a sequence", all_headers, index=0)
        preset_choice = st.selectbox("Quality preset", list(quality_presets.keys()))
        count_choice = st.number_input("Read count", min_value=1, value=100)

        add_btn = st.form_submit_button("Add this sequence")
        if add_btn:
            st.session_state.targets.append({
                "header": seq_choice,
                "number": int(count_choice),
                "sim_quality": preset_choice,
            })
            st.rerun()

    # show selected targets and with option to remove set sequences
    if st.session_state.targets:
        st.write("### Selected Targets")
        for idx, t in enumerate(st.session_state.targets):
            col1, col2 = st.columns([4, 1])
            with col1:
                st.write(f"- {t['header']} | {t['number']} reads | {t['sim_quality']}")
            with col2:
                if st.button("Remove", key=f"remove_{idx}"):
                    st.session_state.targets.pop(idx)
                    st.rerun()

    # set output directory & additional options
    output_dir = st.text_input("Output directory", value="output_LRS")
    run_nanoplot = st.checkbox("Run NanoPlot on combined reads")
    plot_composition = st.checkbox("Plot sample composition")
    remove_temp = st.checkbox("Remove temporary files (simulated sequences)")

    if st.button("Run Simulation") and st.session_state.targets:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, "temporary"), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "logs"), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "temporary", "simulated"), exist_ok=True)

        # extract selected sequences
        found_targets = extract_sequences(fasta_path, st.session_state.targets, output_dir)

        # simulate reads
        for target in found_targets:
            fasta_seq_path = os.path.join(output_dir, "temporary", f"{target['header'].replace(' ', '_')}.fasta")
            cmd = build_command(fasta_seq_path, target["number"], target["sim_quality"], quality_presets)
            output_file = os.path.join(output_dir,"temporary","simulated",
                f"{target['header'].replace(' ', '_')}_{target['sim_quality']}_{target['number']}.fastq")
            
            log_file = os.path.join(output_dir, "logs", f"{target['header'].replace(' ', '_')}_simulation_log.txt")
            print("simulating with command:", cmd)
            simulate_badreads(cmd, output_file, log_file)

        # combine and sample
        combined_output = os.path.join(output_dir, "combined_simulated_reads.fastq")
        combine_and_sample_reads(found_targets, output_dir, combined_output)

        # NanoPlot & composition
        Nanoplot(run_nanoplot, sample_file=combined_output, output_dir=os.path.join(output_dir, "nanoplot"))
        plot_sample_composition_plotly(plot_composition, found_targets, output_dir)

        if remove_temp:
            shutil.rmtree(os.path.join(output_dir, "temporary"), ignore_errors=True)

        st.success("Simulation complete! Check output directory for results.")
