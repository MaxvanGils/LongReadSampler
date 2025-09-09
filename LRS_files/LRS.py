from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import json
import argparse
import subprocess
from random import choices
import plotly.graph_objects as go
import plotly.express as px
import plotly.subplots as sp
import shutil

def load_targets(file_path):
    """
    Load targets from a file with format:
    header <spaces> number <spaces> sim_quality
    """
    targets = []
    with open(file_path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                # split line into parts (3, starting from the right)
                parts = line.strip().rsplit(maxsplit=2)
                header = parts[0]
                number = int(parts[1])
                sim_quality = parts[2]
                targets.append({'header': header, 'number': number, 'sim_quality': sim_quality})
    return targets


def load_quality_presets(file_path):
    """
    Load quality presets used by Badread to simulate reads 
    """
    with open(file_path) as f:
        return json.load(f)

def extract_sequences(fasta_file, targets, output_dir):
    """
    Extract sequences based on headers in targets.
    """
    print("/n" )
    found_targets = set()
    found_target_obj = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        for target in targets:
            if record.description == target['header']:
                found_targets.add(target['header'])
                found_target_obj.append(target)
                name = record.description.replace(" ", "_")
                output_path = os.path.join(output_dir, "temporary", f"{name}.fasta")
                # rearrange the record to have a clean header
                new_record = SeqRecord(
                    seq=record.seq,
                    id=name,             
                    description=name     
                )

                # extract each sequence to its own fasta file
                with open(output_path, "w") as out_handle:
                    SeqIO.write(new_record, out_handle, "fasta")
                print(f"Extracted {record.description} → {output_path} (Number: {target['number']}, Quality: {target['sim_quality']})")

    # Print missing targets
    missing_targets = {t['header'] for t in targets} - found_targets
    for missing in missing_targets:
        print(f"Warning: Target '{missing}' not found in the FASTA file.")
        # write a warning to a log file
        os.makedirs(os.path.join(output_dir, "logs"), exist_ok=True)
        with open(os.path.join(output_dir, "logs", "missing_targets_log.txt"), "a") as log_handle:
            log_handle.write(f"Warning: Target '{missing}' not found in the FASTA file.\n")
    
    return found_target_obj

def build_command(fasta_path, reads, sim_quality, quality_presets):
    """
    Build the badread command based on sim_quality presets
    """
    command = ['badread', 'simulate', '--reference', fasta_path, '--quantity', f"{reads}x"]

    if sim_quality in quality_presets:
        for key, value in quality_presets[sim_quality].items():
            command += [key, value]
    else:
        print(f"Warning: sim_quality '{sim_quality}' not found in presets. Using default command.")
        # add warning to a log file
        with open(os.path.join(output_dir, "logs", "missing_quality_presets_log.txt"), "a") as log_handle:
            log_handle.write(f"Warning: sim_quality '{sim_quality}' not found in presets. Using default command.\n")

    return command

def simulate_badreads(command, output_file, log_file):
    """
    Run the badread command and save output to a file, logging errors.
    """
    with open(output_file, 'w') as output_handle, open(log_file, 'w') as log_handle:
        process = subprocess.Popen(command, stdout=output_handle, stderr=subprocess.PIPE)

        # wait for the process to finish and capture stderr
        stderr_output = process.stderr.read()

        # log any errors to the corresponding log file
        if stderr_output:
            log_handle.write(stderr_output.decode())

    print(f"Generated {output_file} and logged info to {log_file}\n")


def combine_and_sample_reads(found_targets, output_dir, combined_output):
    """
    Combine simulated reads from separate files, sampling to requested count.
    """
    all_records = []
    for target in found_targets:
        # set path to the simulated fastq file
        fastq_path = os.path.join(output_dir, "temporary","simulated", f"{target['header'].replace(' ', '_')}_{target['sim_quality']}_{target['number']}.fastq")
        print(f"Processing {fastq_path} for sampling...")
        requested = target['number'] 
        records = list(SeqIO.parse(fastq_path, "fastq"))
        if not records:
            print(f"Warning: No reads found in {fastq_path}")
            continue

        # if fewer reads are available than requested, take all and sample up to requested
        sampled = choices(records, k=requested)
        all_records.extend(sampled)
        # write how many reads were sampled from each file
        print(f"Sampled {requested} reads from {fastq_path} ({len(records)} available)")
        # write to a log file
        with open(os.path.join(output_dir, "logs", "sampling_log.txt"), "a") as log_handle:
            log_handle.write(f"Sampled {requested} reads from {fastq_path} ({len(records)} available)\n")

    with open(combined_output, "w") as out_handle:
        SeqIO.write(all_records, out_handle, "fastq")
    print(f"Combined sampled reads written to {combined_output}\n")



# as a additional feature, run a NanoPlot on the combined reads
def Nanoplot(run, sample_file, output_dir):
    if run:
        nanoplot_command = (
        f"NanoPlot --fastq {sample_file} --only-report --tsv_stats \
          --info_in_report --N50 --plots kde --outdir {output_dir}")
        print("Performing NanoPlot on the combined reads\n")
        subprocess.run(nanoplot_command, shell=True)

# as an additional feature, plot the composition of the final sample using Plotly
def plot_sample_composition_plotly(plot, found_targets, output_dir):
    if plot:
        # set up data for pie chart and table
        labels = [target['header'] for target in found_targets]
        sizes = [target['number'] for target in found_targets]
        qualities=[target['sim_quality'] for target in found_targets]
        total = sum(sizes)
        percentages = [f"{(n/total)*100:.1f}%" for n in sizes]

        # create subplots: 2 rows, 1 column
        fig = sp.make_subplots(
            rows=2, cols=1,
            specs=[[{'type':'domain'}], [{'type':'table'}]],
            #subplot_titles=["Sample Composition Pie Chart", "Sample Composition Table"],
            row_heights=[0.6, 0.4],
            vertical_spacing=0.06
        )

        # logic for pie chart
        fig.add_trace(
            go.Pie(labels=labels, values=sizes, hole=0.3, textinfo='percent'),
            row=1, col=1
        )

        # logic for the table
        fig.add_trace(
            go.Table(
                header=dict(values=["Species", "Frequency", "Percentage", "Quality"], 
                            font=dict(size=14),
                            height=30),
                cells=dict(values=[labels, sizes, percentages, qualities],
                           font=dict(size=12),
                           height=25),
                columnwidth=[200, 100, 100, 100],
            ),
            row=2, col=1
        )
        
        # update layout for better appearance
        fig.update_layout(
            title=dict(text="Sampled Reads Composition",
                x=0.5,             
                xanchor="center",  
                font=dict(size=20) 
                ),
            height=800,
            template="plotly_white"
        )

        chart_path = os.path.join(output_dir, "sample_composition.html")
        fig.write_html(chart_path)
        print(f"Interactive pie chart and table saved to {chart_path}\n")

if __name__ == "__main__":
    #set up arguments to parse
    argparse = argparse.ArgumentParser(description="Simulate sequence data using Badread and sample data into one file.")
    argparse.add_argument("--fasta", "--F", required=True, help="Input FASTA file with sequences.")
    argparse.add_argument("--targets", "--T", required=True, help="File with target headers, read numbers to simulate and quality presets from the JSON file.")
    argparse.add_argument("--quality", "--Q", default="quality_presets.json", help="JSON file with quality presets.")
    argparse.add_argument("--output", "--O", required=True, help="Output directory for extracted FASTA files.")
    argparse.add_argument("--nanoplot", "--N", action='store_true', help="Run NanoPlot on the combined simulated reads.")
    argparse.add_argument("--plot", "--P", action='store_true', help="Generate a pie chart of the sampled reads composition.")
    argparse.add_argument("--clenanup", "--C", choices=["none", "temp"], default="none", help="Remove intermediate files: 'none' (default), 'temp' (remove temporary dir with simulated sequences).")
    args = argparse.parse_args()
    
    #parse arguments
    fasta_file = args.fasta
    target_file = args.targets
    quality_file = args.quality
    output_dir = args.output
    run = args.nanoplot
    plot = args.plot
    cleanup = args.clenanup
    
    # make output dirs
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(os.path.join(output_dir, "temporary")):
        os.makedirs(os.path.join(output_dir, "temporary"))
    # load targets ans quality presets
    targets = load_targets(target_file)
    quality_presets = load_quality_presets(quality_file)
    # check for which sequences are present in the fasta file and report missing
    found_targets = extract_sequences(fasta_file, targets, output_dir)

    # build Badread command
    print('\n')
    for target in found_targets:
        fasta_path = os.path.join(output_dir, "temporary", f"{target['header'].replace(' ', '_')}.fasta")
        reads = target['number']
        sim_quality = target['sim_quality']
        cmd = build_command(fasta_path, reads, sim_quality, quality_presets)
        print(f"Simulating sequence data for {target} with the following command:", cmd)

        # generate file names for output and log of simulation
        base_name = f"{target['header'].replace(' ', '_')}_{sim_quality}_{reads}"
        output_file = os.path.join(output_dir, "temporary", "simulated", f"{base_name}.fastq")
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "logs"), exist_ok=True)
        log_file = os.path.join(output_dir, "logs", f"{base_name}_simulation_log.txt")

        simulate_badreads(cmd, output_file, log_file)
    
    # combine and sample reads after simulation
    combined_output = os.path.join(output_dir, "combined_simulated_reads.fastq")
    combine_and_sample_reads(found_targets, output_dir, combined_output)

    # optionally run NanoPlot on the combined reads
    Nanoplot(run, sample_file=combined_output, output_dir=os.path.join(output_dir, "nanoplot"))

    # optionally plot the composition of sampled reads
    plot_sample_composition_plotly(plot, found_targets, output_dir)

    # cleanup logic
    if cleanup =="temp":
        shutil.rmtree(os.path.join(output_dir, "temporary"), ignore_errors=True)

