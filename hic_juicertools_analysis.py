import subprocess
import sys

# Define required packages
required_packages = ["pandas", "matplotlib", "upsetplot", "Pillow"]

# Check if required packages are installed
missing_packages = [pkg for pkg in required_packages if not subprocess.run([sys.executable, "-m", "pip", "show", pkg], capture_output=True).returncode == 0]

# Install missing packages
if missing_packages:
    print("Installing missing packages...")
    for pkg in missing_packages:
        subprocess.run([sys.executable, "-m", "pip", "install", pkg], check=True)
    print("Packages installed successfully.")

import argparse  # Module for parsing command-line arguments
import os  # Module for interacting with the operating system
import shutil  # Module for high-level file operations
import multiprocessing  # Module for parallel processing
import pandas as pd  # Library for data manipulation and analysis
import matplotlib.pyplot as plt  # Library for creating plots
import matplotlib.ticker as ticker  # Module for customizing tick locators and formatters
from itertools import combinations  # Function to generate combinations
from upsetplot import from_memberships, plot  # Functions for creating upset plots
import matplotlib.pyplot as plt  # Library for creating plots
from PIL import Image  # Library for working with images

# Function to generate eigenvectors
def generate_eigenvectors(hic_files, output_folder):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    for hic_file in hic_files:
        # Extract file name without extension
        file_name = os.path.splitext(os.path.basename(hic_file))[0]
        # Create a subfolder for each Hi-C file
        folder_name = os.path.join(output_folder, f"{file_name}_AB_compartments")
        os.makedirs(folder_name, exist_ok=True)
        # List of chromosomes including X and Y
        chromosomes = list(range(1, 23)) + ['X', 'Y']
        for chromosome in chromosomes:
            # Define output text file path
            txt_file = os.path.join(folder_name, f"{chromosome}", f"AB_compartments_chr{chromosome}.txt")
            os.makedirs(os.path.join(folder_name, str(chromosome)), exist_ok=True)
            # Command to generate eigenvector
            command = f"java -jar juicer_tools.2.20.00.ac.jar eigenvector VC {hic_file} chr{chromosome} BP 500000 -p {txt_file}"
            # Execute the command
            subprocess.run(command, shell=True)

# Function to generate bedgraphs
def generate_bedgraphs(chr_lengths, output_folder):
    for chr, length in chr_lengths.items():
        for folder in os.listdir(output_folder):
            if os.path.isdir(os.path.join(output_folder, folder)):
                for subfolder in os.listdir(os.path.join(output_folder, folder)):
                    # Define input text file
                    file = os.path.join(output_folder, folder, subfolder, f"AB_compartments_chr{chr}.txt")
                    if os.path.isfile(file):
                        # Generate bedgraph files
                        subprocess.run(
                            f"echo -e 'chr{chr}\\t0\\t{length}' | bedtools makewindows -b stdin -w 500000 | paste - {file} > {os.path.join(output_folder, folder, subfolder, f'chr{chr}.bedgraph')}",
                            shell=True,
                            check=True
                        )
                        subprocess.run(
                            f"awk '{{if ($4 > 0) print $1\"\\t\"$2\"\\t\"$3\"\\t\"\".\"\"\\t\"\".\"\"\\t\"\"+\"; else if ($4 < 0) print $1\"\\t\"$2\"\\t\"$3\"\\t\"\".\"\"\\t\"\".\"\"\\t\"\"-\"}}' {os.path.join(output_folder, folder, subfolder, f'chr{chr}.bedgraph')} | bedtools merge -s -d 10 -c 6 -o distinct | awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4}}' > {os.path.join(output_folder, folder, subfolder, f'AB_compartments_chr{chr}.bedpe')}",
                            shell=True,
                            check=True
                        )

# Function to combine bedgraph files
def combine_bedgraph_files(directory, combined_file):
    with open(combined_file, "ab") as outfile:
        for root, dirs, files in os.walk(directory):
            chromosome = os.path.basename(root)
            for filename in files:
                if filename.endswith(".bedgraph") and filename.startswith(f"chr{chromosome}"):
                    filepath = os.path.join(root, filename)
                    with open(filepath, "rb") as infile:
                        shutil.copyfileobj(infile, outfile)

# Function to generate profile plots
def generate_profile_plots(output_folder):
    for folder in os.listdir(output_folder):
        if os.path.isdir(os.path.join(output_folder, folder)):
            for subfolder in os.listdir(os.path.join(output_folder, folder)):
                for bedgraph in os.listdir(os.path.join(output_folder, folder, subfolder)):
                    if bedgraph.endswith(".bedgraph"):
                        # Read bedgraph file into DataFrame
                        df = pd.read_csv(os.path.join(output_folder, folder, subfolder, bedgraph), sep='\t', header=None, names=['chr', 'start', 'end', 'eigenvector'])
                        df.sort_values(by='start', inplace=True)
                        # Assign compartment based on eigenvector values
                        df['compartment'] = df['eigenvector'].apply(lambda x: 'a' if x >= 0 else 'b')
                        # Create bar plot
                        fig, ax = plt.subplots(figsize=(15, 5))
                        for idx, row in df.iterrows():
                            color = 'green' if row['compartment'] == 'a' else 'blue'
                            bar_height = row['eigenvector']
                            ax.bar(row['start'], bar_height, width=row['end'] - row['start'], color=color)
                        ax.set_xlabel(f'Chromosome {df["chr"].iloc[0]}')
                        ax.set_ylabel('Eigenvector')
                        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:.0f}'.format(x)))
                        # Save plot as PNG
                        plt.savefig(os.path.join(output_folder, folder, subfolder, f"{os.path.splitext(bedgraph)[0]}_profile_plot.png"), format='png')
                        plt.close()

# Function to generate alluvial plot
def generate_alluvial_plot(abcomp_folder, hic1, hic2):
    # R code for generating alluvial plot
    r_code = f"""
    # R code for alluvial plot
    """
    # Write R code to a temporary file
    with open("temp_script.R", "w") as r_script:
        r_script.write(r_code)

    # Execute R script
    try:
        subprocess.run(["Rscript", "temp_script.R"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing R script: {e}")

# Function to run HiCCUPS
def run_hiccups(hic_files, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for hic_file in hic_files:
        # Extract file name without extension
        file_name = os.path.splitext(os.path.basename(hic_file))[0]
        # Create folder for each Hi-C file
        folder_name = os.path.join(output_folder, file_name)
        os.makedirs(folder_name, exist_ok=True)
        # Run HiCCUPS
        subprocess.run(f"java -jar juicer_tools.2.20.00.ac.jar hiccups --cpu --threads {multiprocessing.cpu_count()} {hic_file} {folder_name}", shell=True)

# Function to run APA
def run_apa(loops_folder, hic_files):
    for hic_file in hic_files:
        # Extract file name without extension
        file_name = os.path.splitext(os.path.basename(hic_file))[0]
        sample_folder = file_name
        sample_folder_path = os.path.join(loops_folder, sample_folder)
        sample_apa_folder = os.path.join(sample_folder_path, "APA")
        os.makedirs(sample_apa_folder, exist_ok=True)
        bedpe_file = os.path.join(sample_folder_path, "merged_loops.bedpe")
        if os.path.isfile(bedpe_file):
            # Run APA
            subprocess.run(f"java -jar juicer_tools.2.20.00.ac.jar apa --threads {multiprocessing.cpu_count()} {hic_file} {bedpe_file} {sample_apa_folder}", shell=True)

# Function to create upset plot
def upsetplot(loop_files, loop_names,output_folder):
    # Function to determine if two loops are the same
    def is_same_loop(loop1, loop2):
        # Logic to determine if two loops are the same
        return False

    # Function to count common loops
    def count_common_loops(loops, memberships):
        # Logic to count common loops
        return []

    # Dictionary to store loop data
    loops = {}
    for idx, file in enumerate(loop_files):
        df = pd.read_csv(file, sep='\t', header=None, skiprows=2)
        loops[loop_names[idx]] = [{'chr1': row[0], 'x1': row[1], 'x2': row[2], 
                                'chr2': row[3], 'y1': row[4], 'y2': row[5]} for row in df.values]

    # Generate combinations of loop files
    loop_files_combinations = []
    for r in range(1, len(loop_files) + 1):
        loop_files_combinations.extend(combinations(loop_files, r))

    loop_labels = loop_names
    memberships = [[loop_names[i] for i, file in enumerate(loop_files) if file in comb] for comb in loop_files_combinations]

    memberships.insert(0, [])

    counts = count_common_loops(loops, memberships)

    example_df = from_memberships(memberships, data=counts)

    # Plot upset plot
    plot_obj = plot(example_df)

    plt.xticks(range(len(loop_labels)), loop_labels)

    plt.savefig(os.path.join(output_folder, 'upset_plot.png'))

    plt.show()

# Function to run HiCCUPSDiff
def hiccupsdiff(hic1, hic2, hic1_bedpe, hic2_bedpe, hic_files, output_folder):
    # Run HiCCUPSDiff
    subprocess.run(["java", "-jar", "juicer_tools.2.20.00.ac.jar", "hiccupsdiff", "--cpu", "--threads", str(multiprocessing.cpu_count()), "-k", "VC", hic1, hic2, hic1_bedpe, hic2_bedpe, os.path.join(output_folder, "differential_loops")], check=True)
    
    # Run APA for differential loops
    subprocess.run(f"java -jar juicer_tools.2.20.00.ac.jar apa --threads {multiprocessing.cpu_count()} {hic1} HiC_Analysis/Loops/differential_loops/differential_loops1.bedpe HiC_Analysis/Loops/differential_loops/hic1_hic1_apa", shell=True, check=True)
    subprocess.run(f"java -jar juicer_tools.2.20.00.ac.jar apa --threads {multiprocessing.cpu_count()} {hic1} HiC_Analysis/Loops/differential_loops/differential_loops2.bedpe HiC_Analysis/Loops/differential_loops/hic1_hic2_apa", shell=True, check=True)
    subprocess.run(f"java -jar juicer_tools.2.20.00.ac.jar apa --threads {multiprocessing.cpu_count()} {hic2} HiC_Analysis/Loops/differential_loops/differential_loops2.bedpe HiC_Analysis/Loops/differential_loops/hic2_hic2_apa", shell=True, check=True)
    subprocess.run(f"java -jar juicer_tools.2.20.00.ac.jar apa --threads {multiprocessing.cpu_count()} {hic2} HiC_Analysis/Loops/differential_loops/differential_loops1.bedpe HiC_Analysis/Loops/differential_loops/hic2_hic1_apa", shell=True, check=True)
    
    # Generate plots for APA results
    for resolution in ['10000', '25000', '5000']:
        hic1_hic1 = Image.open(f"HiC_Analysis/Loops/differential_loops/hic1_hic1_apa/{resolution}/gw/APA.png")
        hic1_hic2 = Image.open(f"HiC_Analysis/Loops/differential_loops/hic1_hic2_apa/{resolution}/gw/APA.png")
        hic2_hic1 = Image.open(f"HiC_Analysis/Loops/differential_loops/hic2_hic1_apa/{resolution}/gw/APA.png")
        hic2_hic2 = Image.open(f"HiC_Analysis/Loops/differential_loops/hic2_hic2_apa/{resolution}/gw/APA.png")

        fig, axs = plt.subplots(2, 2)

        axs[0, 0].imshow(hic1_hic1)
        axs[0, 0].set_title(f'{os.path.splitext(os.path.basename(hic_files[0]))[0]} HiC, {os.path.splitext(os.path.basename(hic_files[0]))[0]} Loops')
        axs[0, 1].imshow(hic1_hic2)
        axs[0, 1].set_title(f'{os.path.splitext(os.path.basename(hic_files[0]))[0]} Hi-C, {os.path.splitext(os.path.basename(hic_files[0]))[1]} Loops')
        axs[1, 0].imshow(hic2_hic1)
        axs[1, 0].set_title(f'{os.path.splitext(os.path.basename(hic_files[0]))[1]} Hi-C, {os.path.splitext(os.path.basename(hic_files[0]))[0]} Loops')
        axs[1, 1].imshow(hic2_hic2)
        axs[1, 1].set_title(f'{os.path.splitext(os.path.basename(hic_files[0]))[1]} Hi-C, {os.path.splitext(os.path.basename(hic_files[0]))[1]} Loops')

        for ax in axs.flat:
            ax.axis('off')

        fig.text(0.5, 0.04, '{name of hic1}', ha='center')
        fig.text(0.04, 0.5, '{name of hic2}', va='center', rotation='vertical')

        plt.tight_layout()

        plt.savefig(f"HiC_Analysis/Loops/differential_loops/diff_apa_{resolution}.png")

        plt.show()

# Function to generate TADs
def TADs(hic_files):
    output_folder = 'HiC_Analysis/TADs'
    os.makedirs(output_folder, exist_ok=True)
    
    # Iterate over Hi-C files
    for hic_file in hic_files:
        base_name = os.path.splitext(os.path.basename(hic_file))[0]

        # Run arrowhead to identify TADs
        subprocess.run(['java', '-jar', 'juicer_tools.2.20.00.ac.jar', 'arrowhead', '--threads', str(multiprocessing.cpu_count()), hic_file, f'{output_folder}/{base_name}_tads'])  

# Main function to orchestrate the analysis pipeline
def main(hic_files, binsize):
    abcomp_output_folder = "HiC_Analysis/ABcompartments"
    print("generating eigenvectors")
    generate_eigenvectors(hic_files, abcomp_output_folder)
    
    chr_lengths = {
        1: 248956422,
        2: 242193529,
        3: 198295559,
        4: 190214555,
        5: 181538259,
        6: 170805979,
        7: 159345973,
        8: 145138636,
        9: 138394717,
        10: 133797422,
        11: 135086622,
        12: 133275309,
        13: 114364328,
        14: 107043718,
        15: 101991189,
        16: 90338345,
        17: 83257441,
        18: 80373285,
        19: 58617616,
        20: 64444167,
        21: 46709983,
        22: 50818468,
        'X': 156040895,
        'Y': 57227415
    }
    print("generating bedgraphs")
    generate_bedgraphs(chr_lengths, abcomp_output_folder)
    print("generating profile plots")
    generate_profile_plots(abcomp_output_folder)
    
    if len(hic_files) == 2:
        print("generating alluvial plots")
        combined_bedgraph_files = {}

        bedgraph_files = []
        for sample_folder in os.listdir(abcomp_output_folder):
            if os.path.isdir(os.path.join(abcomp_output_folder, sample_folder)):
                combined_file = os.path.join(abcomp_output_folder, sample_folder + ".bedgraph")
                print(combined_file)
                print(os.path.join(abcomp_output_folder, sample_folder))
                bedgraph_files.append(os.path.join(abcomp_output_folder, sample_folder))
                combine_bedgraph_files(os.path.join(abcomp_output_folder, sample_folder), combined_file)
                combined_bedgraph_files[sample_folder] = combined_file

        if len(bedgraph_files) == 2:
            print("yes")
            bedgraph_file1, bedgraph_file2 = bedgraph_files
            generate_alluvial_plot(abcomp_output_folder, bedgraph_file1, bedgraph_file2)

    loops_output_folder = "HiC_Analysis/Loops"
    print("generating loops")
    run_hiccups(hic_files, loops_output_folder)
    
    print("running apa")
    run_apa(loops_output_folder, hic_files)

    hic_folders = [folder for folder in os.listdir(loops_output_folder) if os.path.isdir(os.path.join(loops_output_folder, folder))]
    
    loop_files = []
    loop_names = []
    for sample_folder in os.listdir(loops_output_folder):
        bedpe_file = os.path.join(loops_output_folder, sample_folder, "merged_loops.bedpe")
        if os.path.isfile(bedpe_file):
            loop_files.append(bedpe_file)
            loop_names.append(os.path.splitext(sample_folder)[0])
    
    if len(bedgraph_files) == 2:
        print("generating upsetplot")
        upsetplot(loop_files, loop_names,loops_output_folder)
    
    if len(bedgraph_files) == 2:
        print("generating differential_loops")
        filename_without_extension_0 = os.path.splitext(os.path.basename(hic_files[0]))[0]
        filename_without_extension_1 = os.path.splitext(os.path.basename(hic_files[1]))[0]
        hiccupsdiff(hic_files[0], hic_files[1], f"HiC_Analysis/Loops/{filename_without_extension_0}/merged_loops.bedpe", f"HiC_Analysis/Loops/{filename_without_extension_1}/merged_loops.bedpe", hic_files, "HiC_Analysis/Loops")
    
    print("generating TADs")
    TADs(hic_files)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Hi-C analysis pipeline')
    parser.add_argument('hic_files', nargs='+', help='List of .hic files')
    parser.add_argument('--binsize', type=int, default=500000, help='Binsize for analysis (default: 500000)')
    args = parser.parse_args()
    main(args.hic_files, args.binsize)
