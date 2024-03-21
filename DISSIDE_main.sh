#!/bin/bash
#############################################################
#DISSIDE: Dynamic In Silico Sample Identification for Discrete Evaluation
#Created by: Erick LeBrun, Chien-Chi Lo, Paul Li, and Jessica Salguero
#Contact: elebrun@lanl.gov
#Version: 1.5
#Version Date: 8/9/23
#LANL O4711
#############################################################

# Usage Info
show_info() {
cat << EOF

DISSIDE: Dynamic In Silico Sample Identification for Discrete Evaluation
Version 1.5
Usage: ./scriptname.sh [options]

Info: DISSIDE uses unsupervised learning to select samples that represent an a priori discrete pattern in the data. It further removes "noisy" samples that do not fit discrete patterns well or represent discrete groups below a defined n value. It then estimates the fit and strength of the a priori discrete pattern using both unconstrained and constrained methods.

Options:
	-h		Display this help page.
*	-i	PATH	Input Sample x Feature type data file in csv format. Samples should be columns and features (e.g OTU/ASV/other) should be in rows.
	-dm	Use this flag in addition to -i if input is a csv format square distance matrix. The -d parameter is irrelevant if -dm flag is used. Raw and Cleaned data output files will be distance matrices.
	-am Use this flag in addition to -i if input is a csv format square adjacency matrix. Raw and Cleaned data output files will be adjacency matrices.
	-o	STR	Name for output directory and files (Default = "DISSIDE_out"). A new folder will be created in the current directory.
	-d	STR	Distance metric to use (Default = "bray" for Bray-Curtis for raw data and Default = "jaccard" for Jaccard for Adjacency). Also available all other available in R package vegan.
	-l	STR Linkage method (Default="single"). Single-linkage method is strongly recommended for optimal DISSIDE function. However, all linkage methods from hclust in R are available.
	-n	INT	Minimal number of samples required to retain an a priori discrete grouping (Default = 2).
	-m	INT Override the default maximum k groups to investigate (Default = 0.5 of samples have become singltons). Large sample numbers or complex data can increase computational time non-linearly and investigating fewer k can improve runtime.
	-p	INT Number of permutations to use in testing (Default = 1000).
	-t	INT Number of threads to use for parallelization (Default = Max cores available).

* Indicates a required field

Dependencies:
DISSIDE_ul_R.R (Included)
R base
pandoc
plotly-orca
R packages:
vegan
ggplot2
plotly
htmlwidgets
RColorBrewer
python-kaleido

EOF
}

# Clear variables and set defaults
unset $inputfile
projname="DISSIDE_out"
distmetric="bray"
ngroup=2
perms=1000
maxgroups="NULL"
linkage="single"
startdir=$(pwd)
usedist="raw"
numcores=$(getconf _NPROCESSORS_ONLN)
DIR=$(dirname $0)

# Set environmental variable to remove erroneous processx verbosity.
export PROCESSX_NOTIFY_OLD_SIGCHLD="TRUE"

# Parse Options
while :; do
	case $1 in
		-h)
			show_info
			exit
			;;
		-i)
			inputfile="$2"
			shift;
			;;
		-dm)
			usedist="dist"
			;;
		-am)
			usedist="am"
			distmetric="jaccard"
			;;
		-o)
			projname="$2"
			shift;
			;;
		-d)
			distmetric="$2"
			shift;
			;;
		-l)
			linkage="$2"
			shift;
			;;
		-n)
			ngroup="$2"
			shift;
			;;
		-m)
			maxgroups="$2"
			shift;
			;;
		-p)
			perms="$2"
			shift;
			;;
		-t)
			numcores="$2"
			shift;
			;;
		*)
			break
	esac
	shift
done

# Break if options not fullfilled
if [ -z "$inputfile" ];then
	echo "***WARNING: -i input community file must be provided***"
	echo "See -h for help."
	exit 1
fi

if [ $ngroup -lt 2 ];then
	echo "***WARNING: -n must be >=2 for optimal group selection process and to avoid overfitting***"
	echo "Please increase n."
	exit 1
fi

# Process projname to the startdir and the real projname
startdir="$(dirname $projname)"
projname="$(basename $projname)"
if [[ $startdir != /* ]]
then
	startdir="$PWD/$startdir"
fi

# Create project driectory file structure
echo "Building project file structure: $(date)"
mkdir -p $startdir/$projname
mkdir -p $startdir/$projname/Data_Files
mkdir -p $startdir/$projname/Figures
mkdir -p $startdir/$projname/Figures/Score_Plots
mkdir -p $startdir/$projname/Figures/2D_NMDS_Plots
mkdir -p $startdir/$projname/Figures/3D_NMDS_Plots

# Run DISSIDE
echo "Running DISSIDE: $(date)"
case $usedist in
	raw)
		Rscript $DIR/DISSIDE_ul_R.R $startdir/$projname $inputfile $projname $distmetric $ngroup $maxgroups $perms $numcores $linkage
		;;
	dist)
		Rscript $DIR/DISSIDE_ul_R_dist.R $startdir/$projname $inputfile $projname $distmetric $ngroup $maxgroups $perms $numcores $linkage
		;;
	am)
		Rscript $DIR/DISSIDE_ul_R_adj.R $startdir/$projname $inputfile $projname $distmetric $ngroup $maxgroups $perms $numcores $linkage
		;;
esac

echo "Run complete: $(date)"
echo "Thank you for using DISSIDE!"

exit 0
