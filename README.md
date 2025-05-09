# DBGWAS Network Analysis
Analyze and interpret the De bruijn graphs output of DBGWAS using graph based method from networkx. 

1. Download all the component network in json format.

  `download_json.py <dbgwas_dir> <dest_dir> <chrome_loc> <download_loc>`
  
   - dbgwas_dir: path to the directory containing the component specific visualization file. Directory should end with 'visualisations/components/'
   - dest_dir: path to directory where all the output json files will be saved
   - chrome_loc: path to chrome driver. You can download the chrome drivers [here](https://chromedriver.chromium.org/downloads). Be sure, to match the driver version with the installed browser version.
   - download_loc: download destination of chrome profile. Default chrome setting is /home/usr/downloads/

2. Generate fasta and metadata file from the components.

# Cite

Poudel Saugat, Hyun Jason, Hefner Ying, Monk Jon, Nizet Victor, Palsson Bernhard O (2023) Interpreting roles of mutations associated with the emergence of S. aureus USA300 strains using transcriptional regulatory network reconstruction eLife 12:RP90668

https://doi.org/10.7554/eLife.90668.2
