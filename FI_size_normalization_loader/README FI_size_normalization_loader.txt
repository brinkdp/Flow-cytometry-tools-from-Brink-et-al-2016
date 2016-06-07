FI_size_normalization_loader.m 
is a custom initiation script for the Knijnenburg morphology correction model for flow cytometry data.
This script enables automated high-throughput data analysis of .fcs files by calling said correction model,
extracting the normalized data and performing calculations on the geometrical mean of the histograms on the
FL1-H channel (can be modified to assess also other channels).

#################################################################################
Copyright (C) 2015-2016 Applied Microbiology, Lund University, Sweden

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#################################################################################


# Background
(Cited from S1 Supplementary Information in Brink et al. 2016)

Flow cytometry, not unlike genome sequencing and RNA-seq, generates Big Data due to its single cell analysis over multiple channels. This calls for novel in silico tools for facilitating post-data processing, especially in experimental designs where the number of runs, replicates and conditions start to increase. We have designed Matlab (Release R2015a, The MathWorks, Inc., Natick, MA, US) and Python (v3, The Python Software Foundation, US)  scripts that within minutes processes raw .fcs-files (here outputted by the BD Accuri flow cytometer; Becton-Dickinson, NJ, US), normalizes the fluorescence intensity (FI) on channel FL1-H to cell size by calling the Knijnenburg morphology normalization model, and plots the geometrical mean of each normalized histogram versus time as either scatter or bar plots. Thus we have generated a high-throughput pipeline from data acquisition to post-processing that significantly decreases the hands-on analysis time. Since there is no set limitation in the number of strains, replicates and conditions that can be assess with these custom scripts, they are applicable both to batch and microtiter plate cultivation data. 


# Requirements: 

1. MATLAB

2. Knijnenburg morphology correction model:                              
Knijnenburg, T.A., Roda, O., Wan, Y., Nolan, G. P., Aitchison, J.D., & Shmulevich, I. (2011). A regression model approach to enable cell morphology correction in high-throughput flow cytometry.Molecular systems biology, 7(1).                  


# Installation instructions:

The main script that needs to be called in order to run the pipe line is the FI_size_normalization_loader.m (Matlab). Users will have to download the Knijnenburg model separately according to the author’s instructions and store FI_size_normalization_loader.m in the root of the Knijnenburg model folder.  Our custom Matlab script recognises .fcs-files stored in the following folder hierarchy: /<strain>/<biological replicate>/<conditions>. Filenames must include the sample time (in hours) and the folder names are used to keep track of the files during the processing. In order to facilitate the control of filenames and folder hierarchy, a custom Python script (v3, The Python Software Foundation, US) was also written: batch_name_change_accuri_fc_data.py.


# Running the script
In the Matlab terminal, type:

FI_size_normalization(Data, styler, err, norm_plot)

The input parameters are as follows:

Data = 		The filepath to your experiment data. Example:                      
 		'\Data\microplate' to get to the microplate data                        
 		The data is required to be in Folders according to                      
 		the following structure:                                                
 		 \Strains\Replicates\Conditions\.fcs files       
                        
styler = 	Choose the plotting style; barplots: 'bar',                        
  		or a scatter plot FI vs time: 'time'.                                   

err = 		Whether or not you would like to work with                                       
		the replicates. If err = 1 you will get two outputs:                    
		all_reps_geo, which is the same as with err = 0, and                    
		reps_struct that has all the the means for the replicates.                      
		If err = 1 the file will also produce one                                    
		figure per strain, in the mode you've chosen ('bar' or 'time'),        
		showing the means of the replicates as well as error bars               
		(=+/- standard deviation of the replicates)

norm_plot= 	Gives plots of the raw data histograms versus 
		the normalized histograms when norm_plot=1 




# Example


Extract .fcs data from your flow cytometer
Make a directory tree: /<strain>/<biological replicate>/<conditions>
E.g.	 /s288c/biological_replicate1/glucose20g-l

The script is run by calling the directory /<strain>
You can have multiple subfolders in the <biological replicate> and <conditions> levels

E.g.	 \s288c\biological_replicate1\glucose20g-l	(s288c is the model strain for Saccharomyces cerevisiae)
	 \s288c\biological_replicate1\glucose100g-l
	 \s288c\biological_replicate2\glucose20g-l
	 \s288c\biological_replicate2\glucose100g-l


NB! Do not store any files outside of the subfolders on the <biological replicate> and <conditions> levels. The script is written to count the number of files in a directory with the assumption that all counted files are directories

The .fcs files are stored on the bottom level. It is required that the file names includes the time of sampling, as this parameter will be parsed by the script and used for data storage and plots.

E.g. 	sample1 0h	(the h is required! do not insert whitespace between digit and the h: 0 h will not be recognized in the current version)
	sample1 2h
	sample1 4h
OR	0h		(this is the shortest possible names for the .fcs files)
	2h
	4h	


NB! It is reccomended to improve readability of the .fcs files by annotating also strain name, condition and biological replicate in the file name, but this is not a requirement for the functionality of the script.


Let's assume that you have stored your .fcs files in \Flow_Data\s288c\biological_replicate1\glucose20g-l

In the Matlab terminal:
>>FI_size_normalization('\Flow_Data\s288c', 'time', 1, 0)

And the output will be a scatter plot over time with errorbars (standard deviation between replicates) and no histogram normalization plots.

The output data can also be saved to the Matlab workspace by
>>[all_reps_geo,rep_struct] =FI_size_normalization('\Flow_Data\s288c', 'time', 1, 0)
