/*
 * 1. segment nuclei
 * 2. segment cell
 * 3. select nucleus
 * 
 * radius (r2) for 2x area of r1: r2 = sqrt(2)*r1
*/

// figure out how many dilation result in 2x the area
selectWindow("nucleus_dilation_test.tif");
rename("nucleus");
 
for (i = 0; i < 10; i++) {
  selectWindow("nucleus");
  run("Duplicate...", "title=nucleus-" + i);  
  
  for (j = 1; j <= i; j++) {
    run("Dilate"); 
  }
  
  run("Create Selection");
  run("Measure");
  
  close("nucleus-" + i);
}

v0 = getResult('Area', 0);
for (i = 0; i < nResults(); i++) {
    v = getResult('Area', i);
    setResult('factor', i, v/v0);
}
updateResults();

// result: five dilations should give me twice the original area


// create donut mask
selectWindow("nucleus");
getDimensions(width, height, channels, slices, frames);
newImage("isolated_nucleus", "8-bit black", width, height, 1);
run("Set Scale...", "distance=1.6151 known=1 unit=micron");
 
selectWindow("nucleus");
run("Create Selection");
selectWindow("isolated_nucleus");
run("Restore Selection");
setColor(255);
fill();
run("Select None");

for (i = 1; i <= 5; i++) run("Dilate");

run("Restore Selection");
setColor(0);
fill();
run("Select None");

// fit donut on actin
Stack.setChannel(2);
run("Duplicate...", "title=actin");
setAutoThreshold("MinError dark");
setOption("BlackBackground", true);
run("Convert to Mask");
for (i = 1; i <= 3; i++) run("Erode");
imageCalculator("AND create", "isolated_nucleus","actin");
rename("donut_mask_adjusted");
run("Create Selection");

// measure YAP
selectWindow("nucleus_test_stack.tif");
Stack.setChannel(3);
run("Restore Selection");
run("Measure");

n = roiManager('count');
for (i = 0; i < n; i++) {
    roiManager('select', i);
    // process roi here
}





