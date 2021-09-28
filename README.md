ForkCall - a program to call replication fork locations
=======================================================

ForkCall is a script to analyse TrAEL-Seq datasets to call the locations of replication forks within the data.  It does this by using the relative abundance of reads on the forward and reverse strands to calculate a directionality value, and then looks at regions where this value shows a constant tradjectory, either up or down, to call replication activation or termination regions.

Installation
------------

ForkCall is written in R.  To install it you just need to download the latest release from the git repository and unzip it onto your machine.

The script uses the ```tidyverse``` packages for manipulation and plotting so if you don't have these installed already you will need to add them using:

```install.packages("tidyverse")```

Please note that the code uses the updated options introduced in ```readR``` v2.0 so if you have an older install of readR/tidyverse then you'll likely need to update readR to get things to work.  After you've loaded tidyverse you can check the version of readr using ```sessionInfo()```.  You need to see ```readr_2.0.0``` (or higher).

To update readR just run:

```install.packages("readr")```

Input data
----------

The input data to the program are tab-delimited files containing pre-quantitated data for multiple samples, with a separate file for the counts on the forward and reverse strands.  The standard format uses the structure of a SeqMonk report file, but as long as you have a file with the same columns in it you can use any program to create the input data.

The quantitated values should be globally normalised linear read counts (eg RPM).  The forward and reverse input files should use the same measurement windows.

You can look at ```ExampleData/example_fwd.txt``` for an example of what the data should look like.


Running the example data
------------------------

The script is set up to run the example data.  To run it you need to open a shell, move to the forkcall directory and run:

```
Rscript forkcall.R
```

This will generate results in the ExampleData directory.  Alternative you can open the ```forkcall.R``` script in a program such as RStudio and run it from there, but be aware that you will need to set the working directory for it to be able to find the example data.


Running your own data
---------------------

To run your own data simply edit the header section of the ```forkcall.R``` script, entering your own filenames and analysis preferences.

The options to look at in the analysis are the degree of smoothing applied to the data, and the minimum size of the output blocks.  You can look at the diagnostic plots which the program generates to see how well it is calling your data.  If it's over-calling then you might need to increase the smoothing or min block sizes, if it's missing zones then reducing these values will help.


Output
------

The program generates two types of output.

1. A CSV file of the locations and directions of all of the called zones
2. (Optional) A set of diagnostic plots showing the original data, the smoothed data and the called zones.

![diagnostic plot](https://raw.githubusercontent.com/s-andrews/forkcall/main/ExampleData/DiagPlots/sampleSample2_chr1.png)

Problems?
---------

Please report any problems either directly in the issues section of the github repository (https://github.com/s-andrews/forkcall/issues) or by email to simon.andrews@babraham.ac.uk






