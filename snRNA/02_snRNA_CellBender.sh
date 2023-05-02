# Load module on BioHPC
$ du -sh --/home2/s204365/*/ | sort -rh  #Directories only


$ module purge && module load shared slurm python/3.8.x-anaconda
$ module load gcc && module load gsl
$ module load cuda80/current

$ cd /project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender

# Create a conda environment and then activate it. 
$ conda create -n CellBender python=3.8
$ rm -rf ~/.conda/pkgs/

# List all discoverable environments. 
$ conda info --envs 
$ source activate CellBender


# Install the pytables module:
$ conda install -c anaconda pytables


# PyTables is built on top of the HDF5 library, using the Python language and the NumPy package. It features an object-oriented interface that, combined with C extensions for the performance-critical parts of the code (generated using Python), makes it a fast tool for interactively browsing, processing and searching very large amounts of data. One important feature of PyTables is that it optimizes memory and disk resources so that data takes much less space (specially if on-flight compression is used) than other solutions such as relational object-oriented databases. 

# Install PyTorch. 
(CellBender) $ conda install pytorch torchvision -c pytorch

# Clone this repository and install CellBender:
(CellBender) $ git clone https://github.com/broadinstitute/CellBender.git
(CellBender) $ pip install -e CellBender

# We proceed to run remove-background on the trimmed dataset using the following command(s)

# Sample AO1
$ cellbender remove-background \
--input /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO1/outs/raw_feature_bc_matrix.h5 \
--output /project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO1/AO1_cellbender.h5 \
--expected-cells 20000 \ #Find expected cells by “web_summary.html” from /outs file
--total-droplets-included 50000 \#Choose a number that goes a few thousand barcodes into the “empty droplet plateau”

# Sample AO2
$ cellbender remove-background \
--input /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO2/outs/raw_feature_bc_matrix.h5 \
--output /project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO2/AO2_cellbender.h5 \
--expected-cells 10000 \ #Find expected cells by “web_summary.html” from /outs file
--total-droplets-included 50000 \#Choose a number that goes a few thousand barcodes into the “empty droplet plateau”


# Sample AO3
$ cellbender remove-background \
--input /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO3/outs/raw_feature_bc_matrix.h5 \
--output /project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO3/AO3_cellbender.h5 \
--expected-cells 10000 \ #Find expected cells by “web_summary.html” from /outs file
--total-droplets-included 50000 \#Choose a number that goes a few thousand barcodes into the “empty droplet plateau”


# Sample AO4
$ cellbender remove-background \
--input /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO4/outs/raw_feature_bc_matrix.h5 \
--output /project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO4/AO4_cellbender.h5 \
--expected-cells 10000 \ #Find expected cells by “web_summary.html” from /outs file
--total-droplets-included 50000 \#Choose a number that goes a few thousand barcodes into the “empty droplet plateau”

# Sample AO5
$ cellbender remove-background \
--input /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO5/outs/raw_feature_bc_matrix.h5 \
--output /project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO5/AO5_cellbender.h5 \
--expected-cells 8000 \ #Find expected cells by “web_summary.html” from /outs file
--total-droplets-included 20000 \#Choose a number that goes a few thousand barcodes into the “empty droplet plateau”



# Sample AO6
$ cellbender remove-background \
--input /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO6/outs/raw_feature_bc_matrix.h5 \
--output /project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO6/AO6_cellbender.h5 \
--expected-cells 8000 \ # targeted number of nuclei
--total-droplets-included 20000 \#Choose a number that goes a few thousand barcodes into the “empty droplet plateau”


# Sample AO7
$ cellbender remove-background \
--input /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO4/outs/raw_feature_bc_matrix.h5 \
--output /project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO4/AO4_cellbender.h5 \
--expected-cells 8000 \ # targeted number of nuclei
--total-droplets-included 20000 \#Choose a number that goes a few thousand barcodes into the “empty droplet plateau”


# Sample AO8
$ cellbender remove-background \
--input /project/Neuroinformatics_Core/Konopka_lab/s204365/CellRanger/RNA/AO8/outs/raw_feature_bc_matrix.h5 \
--output /project/Neuroinformatics_Core/Konopka_lab/s204365/CellBender/AO8/AO8_cellbender.h5 \
--expected-cells 8000 \ # targeted number of nuclei
--total-droplets-included 20000 \#Choose a number that goes a few thousand barcodes into the “empty droplet plateau”


####### Note #######
# Again, here we leave out the --cuda flag solely for the purposes of being able to run this command on a CPU. But a GPU is highly recommended for real datasets. 
# The computation will finish within a minute or two (after ~150 epochs). The tool outputs the following files:
# X.h5: An HDF5 file containing a detailed output of the inference procedure, including the normalized abundance  of ambient transcripts, contamination fraction of each droplet, a low-dimensional embedding of the background-corrected gene expression, and the background-corrected counts matrix (in CSC sparse format). 
# X-filtered.h5: Same as above, though, only including droplets with a posterior cell probability exceeding 0.5. 
# X.csv: the list of barcodes with a posterior cell probability exceeding 0.5.
# A PDF summary of the results showing (1) the evolution of the loss function during training, (2) a rank-ordered total UMI plot along with posterior cell probabilities, and (3) a two-dimensional PCA scatter plot of the latent embedding of the expressions in cell-containing droplets. Notice the rapid drop in the cell probability after UMI rank ~500. 

