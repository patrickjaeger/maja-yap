#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".nd2") suffix

//run("Bio-Formats Macro Extensions");
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
  src_img = getTitle();
  zProjectionPlus();
  
  selectWindow(src_img);
  cyto_yap = getCytoplasmicYAP();
  
  selectWindow(src_img);
  identifyNuclei();
  
  selectWindow(src_img);
  getNuclearYAP(out_dirs[0], cyto_yap);
  
  selectWindow(src_img);
  saveOverlay(out_dirs[1]);
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

function selectPattern(pattern) {
  // Select window that contains the pattern anywhere 
  // in the title; no * required
  images = getList("image.titles");

  for (i = 0; i < images.length; i++){
    if (matches(images[i], ".*" + pattern + ".*")) selectWindow(images[i]);
  }
}

function zProjectionPlus() { 
  // Run max. Z-projection on all channels
  src_img = getTitle();
  run("Z Project...", "projection=[Max Intensity]");
  close(src_img);
  rename(src_img);
}

function segmentCytoplasm() {
  src_img = getTitle();

  // actin segmentation  
  selectWindow(src_img);
  Stack.setChannel(2);
  run("Duplicate...", "title=actin");
  setAutoThreshold("MinError dark");
  setOption("BlackBackground", true);
  run("Convert to Mask");
  run("Erode");
  run("Erode");

  // nucleus segmentation
  selectWindow(src_img);
  Stack.setChannel(1);
  run("Duplicate...", "title=nuclei");

  setAutoThreshold("Otsu dark");
  setOption("BlackBackground", true);
  run("Convert to Mask");
  run("Dilate");
  run("Dilate");
  run("Erode");
  run("Erode");
  run("Erode");
  run("Dilate");

  // create cytoplasm mask
  imageCalculator("Subtract create", "actin","nuclei");
  rename("cytoplasm");
  close("actin");
  close("nuclei");
  selectWindow("cytoplasm");
}

function getCytoplasmicYAP() { 
// function description
  src_img = getTitle();
  segmentCytoplasm();
  run("Create Selection");
  
  // measure cytoplasmic yap
  selectWindow(src_img);
  Stack.setChannel(3);
  run("Restore Selection");

  run("Set Measurements...", "area mean integrated median redirect=None decimal=2");
  run("Measure");
  run("Select None");
  close("cytoplasm");

  // extract yap value
  cyto_yap = getResult("Mean", 0);
  run("Clear Results");
  return cyto_yap;
}

function identifyNuclei() {
  // Segment nuclei and add them to the ROI-manager 
  Stack.setChannel(1);
  run("Duplicate...", "title=nuclei");

  setAutoThreshold("Otsu dark");
  setOption("BlackBackground", true);
  run("Convert to Mask");

  run("Dilate");
  run("Dilate");
  run("Erode");
  run("Erode");
  run("Erode");
  run("Dilate");
  
  run("Analyze Particles...", "size=50-Infinity add");
  
  close("nuclei");
}

function getNuclearYAP(out_dir, cytoplasmicYAP) {
  // Measure and save YAP intensity in the nucleus
  run("Set Measurements...", "area mean integrated median redirect=None decimal=2");
  Stack.setChannel(3);
  img_title = getTitle();
  img_title = substring(img_title, 0, img_title.length-4) + "_res.csv";
  run("Duplicate...", "title=YAP");
  
  selectWindow("YAP");
  roiManager("multi-measure measure_all");
  close("YAP");
  selectWindow("Results");
  
  for (i = 0; i < nResults(); i++) {
    v = getResult("Mean", i);
    setResult("normalized_mean", i, v/cytoplasmicYAP);
  }
  updateResults();
  
  saveAs("Results", out_dir + File.separator + img_title);
}

function saveOverlay(out_dir) {
  // Save overlay of segmented nuclei on YAP channel
  Stack.setChannel(3);
  img_title = getTitle();
  output_title = substring(img_title, 0, img_title.length-4) + "_overlay.csv";
  run("Duplicate...", "title=YAP");
  run("RGB Color");
  
  // draw cytoplasm
  selectWindow(img_title);
  segmentCytoplasm();
  run("Create Selection");
  selectWindow("YAP");
  run("Restore Selection");
  run("Colors...", "foreground=yellow background=white selection=yellow");
  run("Line Width...", "line=2");
  run("Draw");
  close("cytoplasm");
  
  // draw nuclei
  roiManager("Show All without labels");
  run("Colors...", "foreground=cyan background=white selection=yellow");
  run("Line Width...", "line=2");
  roiManager("Draw");
  roiManager("Delete");

  saveAs("Tiff", out_dir + File.separator + output_title);
  close(img_title);
}

function clearResults() { 
  // Clear results, manager and close all windows
  run("Clear Results");
  //roiManager("Delete");
  close("*");
}