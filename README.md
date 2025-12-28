# LongReadSampler

This is the repository for LongReadSampler (LRS), a wrapper for the Badread long read simulator installable via Conda. LRS can be utilized via the commandline or using a [Streamlit](https://streamlit.io/) based browser interface.

This wrapper utilizes [Badread](https://github.com/rrwick/Badread) to simulate long-read datasets and utilizes the simulated data to create samples with preset compositions. This can be especially useful when wanting to test pipelines or tools classification tools to ensure they work with long-read data of varying qualities.  

⚠️ Note: This project is not affiliated with or endorsed by the Badread development team and is purely a personal project.  

## Features
- Data simulation tool that can be run via either the commandline or with a graphical interface using Streamlit
- Preset qualities for generating long reads,enabling the user to define what errors may or may not be present in the simulated data
- Custom absolute composition control: after simulating the data, LRS will pick up to the defined number of reads per set sequence, giving absolute control of sample composition
- Options to generate a quality report of the combined simulated data using [NanoPlot](https://github.com/wdecoster/NanoPlot), as well as a piechart and table of the final composition for later reference wit [Plotly](https://plotly.com/)

## Licenses

This project is licensed under the MIT license. See the [LICENSE](LICENSE) for full terms.

The Badread source code is not included in this project, and is licensed separately under the GNU General Public License v3.0.
See the [Badread repository](https://github.com/rrwick/Badread) for full source code and license.

## Installation and usage

1. Download this repository or clone it using git:

    ```bash
    git clone https://github.com/MaxvanGils/LongReadSampler
    ```

2. Navigate to the LongReadSampler files directory:

    ```bash
    cd LongReadSampler/LRS_files
    ```

3. Create and activate the Conda environment:

    ```bash
    conda env create -f LRS_conda_env.yml
    conda activate LRS
    ```

4. Install Badread and NanoPlot separately:

    ```bash
    conda install -c bioconda badread=0.4.1 nanoplot
    ```

5. Run LRS using the commandline or browser interface:

    - Commandline:

        ```bash
        python LRS.py --help
        ```

        _note that for the command line version, beside the input fasta file and quality preset, a "target.txt" file specifying the number of reads and quality per sequence header from the fasta file is also required as shown in the example below in a tab-delimited format :_

        ```bash
            sequence header 1	1000	very_bad
            sequence_header_2	500 pretty_good
        ```

    - Browser interface:

        ```bash
        streamlit run LRS_GUI.py
        ```

        Afterwards, open the provided local URL in your browser and upload a fasta file to get started. It loads the default quality presets from the `quality_presets.json` file, but you can also add your own presets (see below). Additionally, users can upload a tsv or csv file with set percentages per species to be simulated. The "columns" of the files are set in the GUI, after which the user can then specify the total number of reads to simulate, as well as the preset quality. (Note, this will make all of the to be simulated reads to one quality preset. IF you want to utilise multiple presets, you could add them in parts with different tsv files)

For creating custom quality presets for the simulation, add either a new .json file or add a new entry to the `quality_presets.json` file. See the [Badread documentation](https://github.com/rrwick/Badread) for more information on the parameters.

## Disclaimer

This project is a personal endeavor and is not affiliated with or endorsed by the Badread development team.  

This application is currently only tested locally, not on a server based setting, and thus its functionality in a server-based or production environment cannot be guaranteed.

## Citation

If you use LRS in your research or projects, please cite this repository as follows:

Max van Gils, LongReadSampler: A wrapper for long read simulation using Badread, GitHub repository <https://github.com/MaxvanGils/LongReadSampler>, [2025]

Please also cite the original Badread paper and/or repository:

> Wick RR. Badread: simulation of error-prone long reads. _Journal of Open Source Software_. 2019;4(36):1316. doi:[10.21105/joss.01316](https://doi.org/10.21105/joss.01316).
>[Badread Github](https://github.com/rrwick/Badread)

## Roadmap

- Add option to define relative compositions from a set maximum amount of reads instead of absolute read numbers