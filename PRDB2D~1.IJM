//User opens each file individually.  No other files should be open.
// Tweaked specfically for 4ch czi files provided on 06182020

setBatchMode('true')
setBatchMode('hide');

//Get the File Name
input = getTitle();

run("Split Channels");

//Stack from DIC Channel is closed
DicChannel = "C1-"+input;
selectWindow(DicChannel);
close();

//Green Channel stack renamed
GreenChannel = "C2-"+input;
selectWindow(GreenChannel);
rename("Green");
green = getTitle();

// Red Channel Stack renamed
RedChannel = "C3-"+input;
selectWindow(RedChannel);
rename("Red");
red = getTitle();

//Stack from farred Channel is closed
FarredChannel = "C4-"+input;
selectWindow(FarredChannel);
close();


// Preprocessing on each channel...
selectWindow(green);
// Corrects for flat-fielding
// See v0.2 of this code for how to remove isolated 'debris spots' if needed...
run("Subtract Background...", "rolling=50 stack");  // Keep all nuclei information
//Normalizes red and green intensities to same range
//Ensures Max stack and Find Maxima prominence is not bias by absolute intensity
resetMinAndMax();
run("8-bit");

selectWindow(red);
run("Subtract Background...", "rolling=50 stack");  
resetMinAndMax();
run("8-bit");

//Create Max of Red and Green Channels to be used for tracking
imageCalculator("Max create stack", red, green);
selectWindow("Result of "+red);
rename("Maxmerge");
merge = getTitle();

//Smoothing intra-nuclear detail must be based carefully on typical nuclear diameter and how close together two nuclei typically get
run("Gaussian Blur...", "sigma=4.0 stack");

rename("preprocessed");
preproc = getTitle(); //the final preprocessed result

	
//A new image is created for the output of the Find Maxima
selectWindow(preproc);
getDimensions(width, height, channels, slices, frames);
slice = frames; // time = frames!  This line migth need to be changed to slices depending on file metadata
output = "output";
newImage(output, "8-bit black", width, height, frames); 


for (i = 1; i < (slice+1); i++) {
	
	selectWindow(preproc);
	setSlice(i);

	//THE SLICE SPECIFIC Processing

	//First round of local maxima
	//better to be over aggressive here and filter at level of tracks later - if it is dropped here, it is gone for good
	run("Find Maxima...", "prominence=3 output=[Single Points]"); // generates new image
	// At this point some cells are marked by multiple local maxima
	rename("max1");
	//Spatial merging and then second round of local max groups spatially close maxima into a single maxima
	//REQUIRES PLUGIN -> https://imagej.net/MorphoLibJ
	// There is a VERY careful balancing here btwn merging multiple localization within a single cell vs incorrectly merging localizations of nearby cells
	run("Morphological Filters", "operation=Dilation element=Disk radius=10"); // radius merges 'nearby' maxima ; generates new image
	rename("max1dilation");
	// Local max on merged particles creates single localizations on the merged particles
	run("Find Maxima...", "prominence=1.0 output=[Single Points]");
	run("Gaussian Blur...", "sigma=1"); //gives 'spots' fiji will like
	rename("max2");
	
	//Pastes slice into output stack
	run("Select All");
	run("Copy");
	selectWindow(output);
	setSlice(i);
	run("Paste"); //the slice

	close("max1");
	close("max1dilation");
	close("max2");
		
}

close(preproc);

//The Red, Green, and Find Maxima Output are merged into a stack
run("Merge Channels...", "c1=["+red+"] c2=["+green+"] c3=["+output+"] create");
run("Channels Tool...");
Stack.setDisplayMode("grayscale");

//Quick fix if the file opens as a Z stack instead of a T stack
getDimensions(width, height, channels, slices, frames);
slice = slices;
if (slice>1){
	run("Properties...", "channels=3 slices=1 frames=slice unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
}

setBatchMode("exit and display");