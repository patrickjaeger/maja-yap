




src_img = getTitle();
run("RGB Color");
rename("mask_validation");
selectWindow(src_img);
identifyNuclei();

n = roiManager('count');
for (i = 0; i < n; i++) {
  selectWindow(src_img);
  createEmptyMask();
  roiManager('select', i);
  createDonutMask();
  selectWindow(src_img);
  createFinalMask();
  selectWindow(src_img);
  measureCytoskeletalYAP();
}
clearEverything();

// FUNCTIONS
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

function measureNuclearYAP(out_dir, cytoplasmicYAP) {
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

function measureCytoskeletalYAP() { 
  Stack.setChannel(3);
  img_title = getTitle();
  img_title = substring(img_title, 0, img_title.length-4) + "_cytoplasmicYAP.csv";
  validation_title = substring(img_title, 0, img_title.length-4) + "_segmentation.tif";
  
  // measure YAP
  run("Restore Selection");
  run("Measure");
  run("Select None");
  
  // draw donut
  selectWindow("mask_validation");
  run("Restore Selection");
  run("Draw", "slice");
  
  close("final_mask");
  
  // save results and validation image
  saveAs("Results", out_dir + File.separator + img_title);
  run("Clear Results");
  
  selectWindow("mask_validation");
  rename(validation_title);
  saveAs("Tiff", out_dir + File.separator + validation_title);
  close(validation_title);
}

function clearEverything() { 
  run("Select None");
  roiManager("Delete");
  close("*");
}

