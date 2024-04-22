# HiC-Pipelines
A collection of R, Python, and command line pipelines designed for analyzing and comparing HiC data. These pipelines facilitate the identification and visualization of AB compartments, loops, TADs.

# Pre-Requirements
- Python
- Pip
- R

# Instructions
In order to execute these pipelines, you will need one to two HiC files. In order to utilize the comparison portion of the pipelines, two HiC files will be necessary. Additionally, you will need the juicer tools jar file.

# hic_juicertools_analysis.py
```
python hic_juicertools_analysis.py <path_to_hic1_file> [path_to_hic2_file] [bin_size]
```
## Important
- The second hic sample is optional
- Bin size is also optional; default is 500,000 bp
## Output
### AB Compartments
- Bedgraph file with AB Compartments for each HiC sample provided
```
HiC_Analysis/ABcompartments/sample_name_AB_compartments.bedgraph
```
- Bedgraph files with AB Compartments for each chromosome for each HiC sample provided
```
HiC_Analysis/ABcompartments/sample_name_AB_compartments/chr#/chr#.bedgraph
```
- AB Compartmentalization profile plots for each chromosome for each HiC sample provided
```
HiC_Analysis/ABcompartments/sample_name_AB_compartments/chr#/chr#_profile_plot.png
```
![image](https://github.com/rujulapradeep/HiC-Pipelines/assets/132700660/ec988afb-62c8-4eeb-8182-46484bf11bd3)

- Alluvial Plot showing differences in AB Compartmentalization if two HiC samples provided
```
HiC_Analysis/ABcompartments/alluvial_plot.png
```
![image](https://github.com/rujulapradeep/HiC-Pipelines/assets/132700660/b9bfe283-9b22-40f3-9376-e77f093f5067)


### Loops
- Bedpe file with loops for each HiC sample provided
```
HiC_Analysis/Loops/sample_name/merged_loops.bedpe
```
- APA for loops for each HiC sample provided at resolutions 5000, 10000, 25000
```
HiC_Analysis/Loops/sample_name/APA/resolution#/gw/APA.png
```
![image](https://github.com/rujulapradeep/HiC-Pipelines/assets/132700660/4498293f-341c-4e28-8def-2921d933066b)
- Upset plot showing differences in loops if two HiC samples provided
```
HiC_Analysis/Loops/upset_plot.png
```
![image](https://github.com/rujulapradeep/HiC-Pipelines/assets/132700660/b2460f96-2e0c-4e4f-9e4a-9663f8d9c57a)

- Bedpe files of loops unique to each HiC sample if two HiC samples provided
```
HiC_Analysis/Loops/differential_loops/differential_loops#.bedpe
```
- APA visualizing differential loops at resolutions 5000, 10000, 25000 when two HiC samples provided
```
HiC_Analysis/Loops/differential_loops/diff_apa_resolution#.png
```
![image](https://github.com/rujulapradeep/HiC-Pipelines/assets/132700660/71cd24f5-5b5e-490d-b350-14e62f311192)

### TADs
- File with TADs for each HiC sample provided
```
HiC_Analysis/TADs/
```
- visualization will be added soon
