#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".nd2") suffix
#@ Integer (label = "Nuclear channel", value=1, min=1, max=3, style="spinner") nuclear_ch
#@ Integer (label = "Measurement channel", value=2, min=1, max=3, style="spinner") measurement_ch
#@ Integer (label = "Donut dilations", value=8, min=1, max=20, style="spinner") donut_dilations

// Create output directories if they do not already exist
resOutDir = input + "/.." + File.separator + "results_csv";
imgOutDir = input + "/.." + File.separator + "results_img";
  
if (!File.exists(input + "/.." + File.separator + "results_csv")) {
  File.makeDirectory(resOutDir);
  File.makeDirectory(imgOutDir);
}

// Main
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
  
  // generate duplicate to draw outlines
  selectWindow(src_img);
  Stack.setChannel(measurement_ch);
  run("Duplicate...", "title=mask_validation");
  run("Enhance Contrast", "saturated=0.35");
  run("RGB Color");
  
  // segment entire cells
  segmentCells(2);
  
  // add nuclei to ROImanager and measure nuclear YAP
  selectWindow(src_img);
  identifyNuclei();

  selectWindow(src_img);
  measureNuclei(resOutDir, measurement_ch);

  n = roiManager('count');
  for (j = 0; j < n; j++) {
    selectWindow(src_img);
    createEmptyMask();
    roiManager('select', j);
    createDonutMask();
    selectWindow(src_img);
    createFinalMask();
    selectWindow(src_img);
    measureCytoplasma();
  }
  
  selectWindow(src_img);
  saveCytoplasmicResults(resOutDir, src_img);
  saveOverlay(imgOutDir, src_img);
  clearEverything();
}


// FILE/DIRECTORY FUNCTIONS -----------------------------------

function BFopen(input, file) { 
// Specific for 3-channel Z-stacks.
// 1. Open image using the bioformats importer and set channels to blue, green, red.
// 2. Do Z-projection and close source file.
  
  run("Bio-Formats Importer", 
  "open=[" + input + File.separator + file + "] autoscale color_mode=Custom " +
  "rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT " + 
  "series_0_channel_0_red=0 series_0_channel_0_green=0 series_0_channel_0_blue=255 " + 
  "series_0_channel_1_red=0 series_0_channel_1_green=255 series_0_channel_1_blue=0 " +
  "series_0_channel_2_red=255 series_0_channel_2_green=0 series_0_channel_2_blue=0");
  
  src_img = getTitle();
  run("Z Project...", "projection=[Max Intensity]");
  close(src_img);
  rename(src_img);
}


// MAIN FUNCTIONS ----------------------------------------------------------------

function identifyNuclei() {
  // Segment nuclei and add them to the ROI-manager 
  Stack.setChannel(nuclear_ch);
  run("Duplicate...", "title=nuclei");
  
  // Threshold
  // setOption("ScaleConversions", true);  // needed?
  run("8-bit");
  run("Median...", "radius=2");
  run("adaptiveThr ", "using=Mean from=55 then=-15 output");
  run("Invert");
  run("Invert LUT");
  run("Close-");
  run("Fill Holes");
  run("Analyze Particles...", "size=40-Infinity circularity=0.3-1.00 add");
  close("nuclei");
  
  // Make new mask without small particles.
  // Advantage: erode/dilate can be performed without 
  //            interference from small particles
  newImage("nuclei2", "8-bit black", 1024, 1024, 1);
  setColor(255);
  n = roiManager('count');
  for (i = 0; i < n; i++) {
    roiManager("Select", i);
    run("Fill", "slice");
  }
  run("Select None");
  roiManager("reset");
  close("Mask");
  
  // Add only large particles to manager
  selectWindow("nuclei2");
  run("Erode");
  run("Analyze Particles...", "size=40-Infinity circularity=0.3-1.00 add");
  close("nuclei2");
  
  // draw outlines of nuclei on RGBduplicate
  selectWindow("mask_validation");
  setForegroundColor(255, 255, 0);
  n = roiManager('count');
  for (i = 0; i < n; i++) {
    roiManager("Select", i);
    run("Draw", "slice");
  }
  roiManager("Deselect");
  
}
//identifyNuclei();

function measureNuclei(out_dir, measurement_channel) {
  // Measure intensity in the nucleus
  Stack.setChannel(measurement_channel);
  run("Set Measurements...", "area mean integrated median redirect=None decimal=2");
  roiManager("Measure");
  
  img_title = getTitle();
  output_title = substring(img_title, 0, img_title.length-4) + "_nuclear-signal.csv";
  //selectWindow("Results");
  saveAs("Results", out_dir + File.separator + output_title);
  run("Clear Results");
}

//// Donut segmentation ------
function segmentCells(channel) { 
// function description
  
  //// threshold
  selectWindow(src_img);
  Stack.setChannel(channel);
  run("Duplicate...", "title=cell_segm");
  run("8-bit");
  run("Median...", "radius=2");
  run("adaptiveThr ", "using=Mean from=101 then=-2 output");
  selectWindow("Mask");
  run("Analyze Particles...", "size=100-Infinity show=Overlay add");
  close("cell_segm");
   
  //// draw large particles into new image
  newImage("cell_mask", "8-bit black", 1024, 1024, 1);
  setColor(255);
  n = roiManager('count');
  for (i = 0; i < n; i++) {
    roiManager("Select", i);
    run("Fill", "slice");
  }
  run("Select None");
  roiManager("reset");
  run("Erode");
  run("Create Selection");
  close("Mask");
  
  // draw outline for validation
  selectWindow("mask_validation");
  run("Restore Selection");

  setForegroundColor(0, 255, 255);  // yellow
  run("Draw", "slice");
}
//segmentCells(2);

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
  for (i = 1; i <= donut_dilations; i++) run("Dilate");
  
  //// cut out donut hole
  run("Restore Selection");
  setColor(0);
  fill();
  run("Select None");
}

function createFinalMask() { 
// isolate overlap between donut and cell to create donut
  Stack.setChannel(measurement_ch);
  imageCalculator("AND", "isolated_nucleus", "cell_mask");
  selectWindow("isolated_nucleus");
  run("Create Selection");
  
}

function measureCytoplasma() {
  
  
  // check that image is not empty
  selectWindow("isolated_nucleus");
  getDimensions(width, height, channels, slices, frames);
  px_sum = 0;
  for (x = 1; x <= width; x++) {
    for (y = 1; y <= height; y++) {
      px_xy = getPixel(x, y);
      px_sum = px_sum + px_xy;
    }
  }
  
  // measure cytoplasma
  selectWindow(src_img);
  Stack.setChannel(measurement_ch);
  if (px_sum > 0) run("Restore Selection");
  run("Measure");
  run("Select None");
  
  // draw donut for validation
  selectWindow("mask_validation");
  if (px_sum > 0) {
    run("Restore Selection");
    setForegroundColor(255, 0, 255);  // magenta
    run("Draw", "slice");
  }
  close("isolated_nucleus");
}

function saveCytoplasmicResults(csv_out_dir, img_title) {
  // Save results and validation image
  results_title = substring(img_title, 0, img_title.length-4) + "_cytoplasma-signal.csv";
  saveAs("Results", csv_out_dir + File.separator + results_title);
  run("Clear Results");
}

function saveOverlay(img_out_dir, img_title) {
  // Save validation image
  validation_title = substring(img_title, 0, img_title.length-4) + "_segmentation.tif"; 
  
  selectWindow("mask_validation");
  saveAs("Tiff", img_out_dir + File.separator + validation_title);
  close(validation_title);
}

function clearEverything() { 
  run("Select None");
  roiManager("reset");
  run("Clear Results");
  close("*");
}
//clearEverything();