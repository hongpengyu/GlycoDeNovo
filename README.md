# GlycoDeNovo in Java

Reference: GlycoDeNovo - an Efficient Algorithm for Accurate de novo Glycan Topology Reconstruction from Tandem Mass Spectra. (https://www.ncbi.nlm.nih.gov/pubmed/28786094)

# Notes:

- CGlycoDeNovo.java is the main part of algorithm implementation in the paper, including interpret peaks and reconstruct
- CSpectrum.java is the part of reading, processing data file and outputing the result of reconstruction
- Use https://github.com/lbehnke/hierarchical-clustering-java for peak single linkage hierarchical clustering with cutoff based on distance.
- Sort the result in peaks by Collection.sort according to its score(size) and Default String order
- Please see more detailed comments in codes


# To Run the project :

Using intellj IDE to run test.java or Using commandline to run the precompiled GlycoDeNovo.jar

#### You should specify the following parameters by yourself
> TEST_DATA_PATH TEST_RES_PATH check2H checkGap PPM
#### Default Setting
>  DATA_PATH = ".\\tests\\data\\";
   RES_PATH = ".\\tests\\results\\";
   check2H = false;
   checkGap = false;
   PPM = 5;
          

#### Run GlycoDeNovo using Command Line

You can specify input with commandline to run the precompile executable: GlycoDeNovo.jar

Run code with default setting :

> java -jar GlycoDeNovo.jar

Run with new settings:

> java -jar GlycoDeNovo.jar DATA_PATH=.\\tests\\data\\A2F.DR.Na.EED.PM.txt RES_PATH=.\\tests\\results\\ checkGap=false check2H=true PPM=5

Run code with specific condition (DATA_PATH RES_PATH check2H checkGap AccuracyPPM)

Note that the DATA_PATH can also be a folder if you want to run with all the data inside that folder
> java -jar GlycoDeNovo.jar DATA_PATH=.\\tests\\data\\ RES_PATH=.\\tests\\results\\

You do not have to specify all the params, you can just specify params you need to change

> java -jar GlycoDeNovo.jar check2H=true PPM=10

> java -jar GlycoDeNovo.jar DATA_PATH=.\\tests\\data\\A2F.DR.Na.EED.PM.txt

> java -jar GlycoDeNovo.jar checkGap=true

If you want to run this code with IntelliJ, you should run with test.java in src folder.
The input command format is similar to run GlycoDeNovo.jar using commandline.
