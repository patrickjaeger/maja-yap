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
  // open image and do Z projection
  BFopen(input, file);
  src_img = getTitle();
  zProjectionPlus();
  
  // generate duplicate to draw outlines
  selectWindow(src_img);
  Stack.setChannel(3);
  run("Duplicate...", "title=mask_validation");
  run("RGB Color");
  
  // add nuclei to ROImanager and measure nuclear YAP
  selectWindow(src_img);
  identifyNuclei();
  
  selectWindow(src_img);
  measureNuclearYAP(out_dirs[0]);

  n = roiManager('count');
  for (i = 0; i < n; i++) {
    selectWindow(src_img);
    createEmptyMask();
    roiManager('select', i);
    createDonutMask();
    selectWindow(src_img);
    createFinalMask();
    selectWindow(src_img);
    measureCytoplasmicYAP();
  }
  
  selectWindow(src_img);
  saveCytoplasmicResults(out_dirs[0], out_dirs[1]);
  clearEverything();
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
  //run("Dilate");
  
  run("Analyze Particles...", "size=50-Infinity add");
  
  close("nuclei");
}

function measureNuclearYAP(out_dir) {
  // Measure and save YAP intensity in the nucleus
  Stack.setChannel(3);
  img_title = getTitle();
  img_title = substring(img_title, 0, img_title.length-4) + "_nuclearYAP.csv";
  run("Duplicate...", "title=YAP");
  
  selectWindow("YAP");
  run("Set Measurements...", "area mean integrated median redirect=None decimal=2");
  roiManager("multi-measure measure_all");
  close("YAP");
  
  selectWindow("Results");
  saveAs("Results", out_dir + File.separator + img_title);
  run("Clear Results");
}

function createEmptyMask() { 
  // create empty image
  getDimensions(width, height, channels, slices, frames);
  getPixelSize(unit, pixelWidth, pixelHeight);
  newImage("isolated_nucleus", "8-bit black", width, height, 1);
  run("Set Scale...", "distance=" + 1/pixelWidth + " known=1 unit=micron");
}

function createDonutMask() { 
  //// draw nucleus in new window and inflate
  selectWindow("isolated_nucleus");
  setColor(255);
  fill();
  run("Select None");
  for (i = 1; i <= 5; i++) run("Dilate");
  
  //// cut out donut hole
  run("Restore Selection");
  setColor(0);
  fill();
  run("Select None");
}

function createFinalMask() { 
  // intersection of donut and actin
  Stack.setChannel(2);
  run("Duplicate...", "title=actin");
  setAutoThreshold("MinError dark");
  setOption("BlackBackground", true);
  run("Convert to Mask");
  for (i = 1; i <= 3; i++) run("Erode");
  imageCalculator("AND create", "isolated_nucleus","actin");
  rename("final_mask");
  run("Create Selection");
  
  close("actin");
  close("isolated_nucleus");
}

function measureCytoplasmicYAP() { 
  Stack.setChannel(3);
  
  // measure YAP
  run("Restore Selection");
  run("Measure");
  run("Select None");
  
  // draw donut
  selectWindow("mask_validation");
  run("Restore Selection");
  run("Draw", "slice");
  
  close("final_mask");
}

function saveCytoplasmicResults(csv_out_dir, img_out_dir) {
  // save results and validation image
  img_title = getTitle();
  img_title = substring(img_title, 0, img_title.length-4) + "_cytoplasmicYAP.csv";
  validation_title = substring(img_title, 0, img_title.length-4) + "_segmentation.tif"; 
  
  saveAs("Results", csv_out_dir + File.separator + img_title);
  run("Clear Results");
  
  selectWindow("mask_validation");
  rename(validation_title);
  saveAs("Tiff", img_out_dir + File.separator + validation_title);
  close(validation_title);
}

function clearEverything() { 
  run("Select None");
  roiManager("Delete");
  close("*");
}