#@ File (style="open") inputFile

// open and create Z-projection
BFopen(inputFile);
src_img = getTitle();
run("Z Project...", "projection=[Max Intensity]");
close(src_img);
rename(src_img);

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

// cytoplasm mask
imageCalculator("Subtract create", "actin","nuclei");
rename("cytoplasm");
close("actin");
close("nuclei");
selectWindow("cytoplasm");
run("Create Selection");
close("cytoplasm");

// measure cyplasmic yap
selectWindow(src_img);
Stack.setChannel(3);
run("Restore Selection");

run("Set Measurements...", "area mean integrated median display redirect=None decimal=1");
run("Measure");

// ...
cyto_yap = getResult("mean", 0);



// FUNCTIONS
function BFopen(input_file) { 
  // Open input using the bioformats importer
  run("Bio-Formats Importer", 
  "open=[" + input_file + 
  "] autoscale color_mode=Default rois_import=[ROI manager]" +
  " view=Hyperstack stack_order=XYCZT");
}

function binarize() { 
// set all non-zero pixels to 255

  getDimensions(width, height, channels, slices, frames);
  for (x = 0; x < width-1; x++) {
   for (y = 0; y < height-1; y++) {
      px = getPixel(x, y);
      //print(px);
      if (px > 0) {
        setPixel(x, y, 255);
      }
   }
  }
  
}