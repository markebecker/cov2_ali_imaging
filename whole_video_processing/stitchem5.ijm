/*
 * Stitch multiple .dv files in a folder
 */
//input = getDirectory("folder1");
//Dialog.create("Stitch folder contents");

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".dv") suffix
#@ Boolean (label = "these puppies flipped the x?", value = checkbox) flipped
#@ Boolean (label = "do you want these projected?", value = checkbox) projected
processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	img = input + File.separator + file;
	print("Stitching: " + input + File.separator + file);
	if(flipped == false) { 
		print("i am running the flipped == false fork");
		run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by image metadata] browse=&input multi_series_file=&img fusion_method=[Linear Blending] regression_threshold=0.05 max/avg_displacement_threshold=0.05 absolute_displacement_threshold=0.05 compute_overlap computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");
	} else {
		print("i am not running the flipped == false fork");
		run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by image metadata] browse=&input multi_series_file=&img fusion_method=[Linear Blending] regression_threshold=0.05 max/avg_displacement_threshold=0.05 absolute_displacement_threshold=0.05 compute_overlap invert_x computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");
		}
	print("it been stitched");
	if(projected == true) {
		run("Z Project...", "projection=[Max Intensity] all");
		print("projected too");
		pathToOutputFile = output + File.separator + file;
		print("Saving to: " + pathToOutputFile);
		saveAs("tiff", pathToOutputFile);
		print("saved it. moving on...");
		run("Close");
		run("Close");
	}
	else {
		pathToOutputFile = output + File.separator + file;
		print("Saving to: " + pathToOutputFile);
		saveAs("tiff", pathToOutputFile);
		print("saved it. moving on...");
		run("Close");
	}
}