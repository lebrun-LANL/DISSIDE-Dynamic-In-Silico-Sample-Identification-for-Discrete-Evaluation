## Changelog (O4711)
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
