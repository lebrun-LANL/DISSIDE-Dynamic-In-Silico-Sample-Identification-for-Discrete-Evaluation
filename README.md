# DISSIDE: Dynamic In Silico Sample Identification for Discrete Evaluation

DISSIDE uses novel unsupervised learning to select samples that best represent the strongest a priori discrete pattern in a given data set. It further removes major outliers and "noisy" samples that do not fit discrete patterns well or represent outliers in groups below a defined n value. It then estimates the fit and strength of the a priori discrete pattern using both unconstrained and constrained methods for raw and cleaned data.

## Software Overview

Dynamic In Silico Sample Identification for Discrete Evaluation (DISSIDE) uses unsupervised machine learning and mathematical algorithms to determine optimal hierarchical clustering results. DISSIDE accepts raw sample x feature data, distance matrices, or adjacency matrices as input. DISSIDE is best run with default parameters and requires no user input beyond the input data but can be heavily parameterized by the user if desired. DISSIDE first identifies samples with no reasonably related samples and removes those samples from the dataset. DISSIDE then confirms that all data is coherent for relational analysis and comparison. If samples exist without coherence for comparison, DISSIDE will automatically partition and output the data into coherent data sets for the user. On coherent data, DISSIDE uses hierarchical clustering to build a tree of sample relationships. DISSIDE forces a range of *k* groupings onto the tree. DISSIDE calculates a unique DISSIDE score for each *k* using PERMANOVA and ANOSIM tests without allowing for overfitting by using a unique algorithm including weighting that allows DISSIDE to know when to stop evaluating *k*'s. The highest DISSIDE score is selected as optimal *k*. Constrained, unconstrained, and MANOVA-type visualization and analysis is performed on the data at optimal *k*. Groups at optimal *k* with fewer than n samples are cleaned from the data as pattern outliers and constrained, unconstrained, and MANOVA-type visualization and analysis is performed on the cleaned data. Centroids/representative samples are also calculated for each cluster during analysis of raw and clean data.

## Included Script Files

DISSIDE\_main.sh	<br />
DISSIDE\_ul\_R.R	<br />
DISSIDE\_ul\_R_adj.R	<br />
DISSIDE\_ul\_R_dist.R	<br />

## Other Included Files

Example test data is included in the folder TEST_DATA <br />
Strong_pattern_PA <br />
> - 9010_otu_no_noise.csv Table represents a strong pattern of presence/absence data with 10 groups with 10 samples in each group. 1000 total unique and fictional OTUs are represented. Each sample has roughly 180 OTUs with 90 OTUs being 90% similar within group and 90 OTUs being 10% similar outside of group.
>
> - 9010_otu_noise_S25.csv Table has the same data as 9010_otu_no_noise.csv but adds 25 completely randomized samples containing 180 OTUs randomly dispersed across the 1000 OTUs.

Strong_pattern_AB <br />
> - Abund_1050_otu_no_noise.csv Table represents a strong pattern of abundance data with 10 groups with 10 samples in each group. 1000 total unique and fictional OTUs are represented. Each sample has 1000 OTUs with random abundance up to 50% variable within group.
>
> - Abund_1050_otu_noise_S25.csv Table has the same data as Abund_1050_otu_no_noise.csv but adds 25 completely randomized samples containing 1000 OTUs with random abundances.

## Installation

Scripts can be run directly from the program directory. All scripts must remain in the same directory.

## Dependencies
[R base](https://www.r-project.org/) <br />
[pandoc](https://pandoc.org/) <br />
[python-kaleido] (https://github.com/plotly/Kaleido) <br />

R packages: <br />
\- [vegan](https://cran.r-project.org/web/packages/vegan/index.html) <br />
\- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) <br />
\- [plotly](https://cran.r-project.org/web/packages/plotly/index.html) <br />
\- [RColorBrewer](https://cran.r-project.org/package=RColorBrewer) <br />
\- [htmlwidgets](https://cran.r-project.org/package=htmlwidgets)
\- [reticulate] (https://cran.r-project.org/package=reticulate) <br /> 

## Usage

```bash
./DISSIDE_main.sh [options]
```
Options: <br />
\-h		Display this help page. <br />
\*	-i	PATH	Input Sample x Feature type data file in csv format. Samples should be columns and features (e.g OTU/ASV/other) should be in rows. <br />
\-dm	Use this flag in addition to -i if input is a csv format square distance matrix. The -d parameter is irrelevant if -dm flag is used. Raw and Cleaned data output files will be distance matrices. <br />
\-am Use this flag in addition to -i if input is a csv format square adjacency matrix. Raw and Cleaned data output files will be adjacency matrices. <br />
\-o	STR	Name for output directory and files (Default = "DISSIDE_out"). A new folder will be created in the current directory. <br />
\-d	STR	Distance metric to use (Default = "bray" for Bray-Curtis for raw data and Default = "jaccard" for Jaccard for Adjacency). Also available all other available in R package vegan. <br />
\-l	STR Linkage method (Default="single"). Single-linkage method is strongly recommended for optimal DISSIDE function. However, all linkage methods from hclust in R are available. <br />
\-n	INT	Minimal number of samples required to retain an a priori discrete grouping (Default = 2). <br />
\-	-m	INT Override the default maximum k groups to investigate (Default = 0.5 of samples have become singltons). Large sample numbers or complex data can increase computational time non-linearly and investigating fewer k can improve runtime. <br />
\-p	INT Number of permutations to use in testing (Default = 1000). <br />
\-t	INT Number of threads to use for parallelization (Default = Max cores available). <br />
<br />
\* Indicates a required field

## Output Files

The software will create a new project directory. <br />
\- Raw data cleaned for distance based outliers and fully cleaned data will be placed in project directory "Data_Files" folder as csv. <br />
\- Distance based hierarchical tree visualization will be placed in the project directory "Figures/Score_Plots" folder as pdf. <br />
\- PERMANOVA and ANOSIM visualization plots that feed into the DISSIDE score will be placed in the project directory "Figures/Score_Plots" folder as pdf and html. <br />
\- DISSIDE scores used to pick optimal *k* will be placed as pdf file in project directory "Figures/Score_Plots" folder as pdg and html. <br />
\- Static and interactive 2D and 3D NMDS visualizations for raw and cleaned data will be placed in the project directory "Figures" folder as pdf and html. <br />
\- Removed sample data will be placed in the project directory "Data_Files" folder as csv. <br />
\- Constrained, unconstrained and MANOVA-type analysis summary file will be placed in the project directory as csv. <br />
\- Individually labelled metadata files that describe the a priori discrete groupings for the raw and cleaned data are placed in the project directory "Data_Files" folder as csv. <br />

## Notes

\- If no outlier or "noisy" samples are detected outside the a priori pattern below n, the software will still provide and analysis on the a priori detected pattern.

## Planned Features

\- DISSIDE CRAN R package.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change. <br />
Please make sure to update tests as appropriate.

## License
This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Contact Information

This software is developed and maintained by the bioinformatics team in B35 Biosecurity and Public Health at Los Alamos National Laboratory. (O4711) <br />
Primary contributors: Erick S Lebrun, PhD; Paul Li; Chien-Chi Lo; Jessica Salguero <br />

Contact: Erick LeBrun elebrun@lanl.gov

## Changelog
\- Version:1.5  Date:8/9/23 <br />
> - Implemented feature for calculating the centroids (representative sample) for each cluster. This feature is utilized during analysis of raw and clean data. 
> - HTMLwidgets still produces an error stating that 'selfcontained is deprecated.' However, conda environment has most up-to-date version and documentation still utilizes selfcontained. Code works and produces images, but the SaveWidget() function calls may need to be updated in the future (newer versions of HTMLwidgets). 

\- Version:1.4  Date:4/27/23 <br />
> - Replaced adonis calls with adonis2. Adonis has become a deprecated function from the vegan package. Fixed resulting errors in code.
> - Replaced orca calls with kaleido, as orca has become deprecated.
> - HTMLwidgets produces an error stating that 'selfcontained is deprecated.' However, conda environment has most up-to-date version and documentation still utilizes selfcontained. Code works and produces images, but the SaveWidget() function calls may need to be updated in the future (newer versions of HTMLwidgets).
> - Unable to update the Dockerfile to work with Kaleido. Will be removing from the package for now. 

\- Version:1.3  Date:3/10/2021 <br />
> - Added flags and functionality support for distance and adjacency matrix inputs.
> - Added 2D NMDS plots for all axis views.
> - Added static and interactive 3D NMDS plots.
> - Cleaned up some redundant code.
> - Introduced some general optimizations around statistical testing and plotting.

\- Version:1.2  Date:3/3/2021 <br />
> - Change variable and file verbiage to Feature instead of OTU as DISSIDE works on any Sample Feature data, not just OTU.
> - Added additional file structure to help organize output files.
> - Added a check that input data is cohesive (has overlapping features).
> - Added a check for obvious initial outliers.
> - If data has samples where distances cannot be accurately calculated, DISSIDE will partition data into new, cohesive data subsets for the user.

\- Version:1.1  Date:2/3/2021 <br />
> - Added a weight element to score calculation. Scores are now contradictorily weighted positively for k both retaining samples as well as excising groups with samples fewer than n. The result is much better optimal k selection.
> - Added optimization where test computation will only be run on k's that result in 2 or more retained groups.
> - Fixed several hard coded variables that should have referenced args that could have resulted in bugs.
> - Added or re-coded various other optimizations and bug fixes.
