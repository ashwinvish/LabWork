
# get current set of layers
layers = Display.getFront().getLayerSet().getLayers()

for layer in layers:
	# Get all images in the current layer
	patches = layer.getDisplayables(Patch)

	# get reference patch
	reference = patches.get( 0 )

	# get transforms
	bounds = reference.getBoundingBox()
  	at = reference.getAffineTransform()

	target = patches.get( patches.size() - 1 )
  	
  	# apply transforms to target
  	target.translate(bounds.x, bounds.y, False)
  	target.setAffineTransform( at )

# repaint display
Display.repaint()
	

