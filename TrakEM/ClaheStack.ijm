getDimensions( width, height, channels, slices, frames );
//Ij.log("Slices= " +slices);
for (i=2;i<=1345;i++)
	{
		Stack.setSlice(i);	
		run("Enhance Local Contrast (CLAHE)", "blocksize=63"
  + " histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
   	  }