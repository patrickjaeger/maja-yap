#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".nd2") suffix

run("Bio-Formats Macro Extensions");
out_dirs = createOutDirs(input);
processFolder(input);

// FUNCTIONS ----

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, list[i]);
	}
}

function processFile(input, file) { 
  // Process file
  BFopen(input, file);
  runProjection();
  addNucsToManager("C1");
  measureYAP(out_dirs[0]);
  saveNuclearOverlay(out_dirs[1]);
  clearResults();
}

function BFopen(input, file) { 
  // Open input using the bioformats importer
  run("Bio-Formats Importer", 
  "open=[" + input + File.separator + file + 
  "] autoscale color_mode=Default rois_import=[ROI manager]" +
  " view=Hyperstack stack_order=XYCZT");
}

function createOutDirs(input) {
  // Create output directories if they do not already exist
  resOutDir = input + "/.." + File.separator + "results_csv";
  imgOutDir = input + "/.." + File.separator + "results_img";
  
  if (!File.exists(input + "/.." + File.separator + "results_csv")) {
    File.makeDirectory(resOutDir);
    File.makeDirectory(imgOutDir);
  }
  
  dirs = newArray(2);
  dirs[0] = resOutDir;
  dirs[1] = imgOutDir;
  
  return dirs;
}
//td = createOutDirs("fu");
//Array.show(td);

function selectPattern(pattern) {
  // Select window that contains the pattern anywhere 
  // in the title; no * required
  images = getList("image.titles");

  for (i = 0; i < images.length; i++){
    if (matches(images[i], ".*" + pattern + ".*")) selectWindow(images[i]);
  }
}

function runProjection() { 
  // Run max. Z-projection on all channels
  source = getList("image.titles");
  run("Split Channels");
  close(source[0]);
  
  images = getList("image.titles");
  for (i = 0; i < images.length; i++) {
    selectWindow(images[i]);
    run("Z Project...", "projection=[Max Intensity]");
    close(images[i]);
  }
}
//runProjection();

function addNucsToManager(pattern) {
  // Segment nuclei and add them to the ROI-manager
  selectPattern(pattern);
  
  setAutoThreshold("Otsu dark");
  setOption("BlackBackground", true);
  run("Convert to Mask");
  
  run("Fill Holes");
  run("Open");
  run("Watershed");
  
  run("Analyze Particles...", "size=50-Infinity add");
  
  close(pattern);  
}
//addNucsToManager();

function saveNuclearOverlay(out_dir) {
  // Save overlay of segmented nuclei on YAP channel
  
  selectPattern("C3");
  run("RGB Color");
  roiManager("Show All without labels");
  run("Colors...", "foreground=yellow background=white selection=yellow");
  run("Line Width...", "line=2");
  roiManager("Draw");

  img_title = getTitle();
  img_title = substring(img_title, 7, img_title.length-4) + "_overlay.tif";
  saveAs("Tiff", out_dir + File.separator + img_title);
  close(img_title);
}

function measureYAP(out_dir) {
  // Measure and save YAP intensity in the nucleus
  run("Set Measurements...", "area mean median redirect=None decimal=2");
  selectPattern("C3");
  img_title = getTitle();
  img_title = substring(img_title, 7, img_title.length-4) + "_res.csv";
  
  roiManager("multi-measure measure_all");
  selectWindow("Results");
  saveAs("Results", out_dir + File.separator + img_title);
}
//measureYAP("C:/Users/jaege/OneDrive - ETHZ/Dokumente/PhD/R projects/201118_yaptaz_colloc/test");

function clearResults() { 
  // Clear results, manager and close all windows
  run("Clear Results");
  roiManager("Delete");
  close("*");
}
//clearResults();