%fileName = '/Users/alex/Documents/Dropbox/em_images_ashwin/alex_compare_scripts/skeletons/ashwin/Int1_5 [treeline] #20373.swc';
fileName = [pwd '/skeletons/ashwin/Int1_5 [treeline] #20373.swc'];

initiationSite{5} = [56 75 91 93 104 105 106 107 108 109 110 111 112 114 115 116 113 ];
SpecialNodeType = []; 
validNodeTypes = [-1:5];
uniformCSAreas = [];
highlightedNodes = 1;
newFigure = true;
colorString = {'g'};
pixelUnits = true;
cellnumber = 5;

tree = generateIrreducibleDoubleLinkedTree(fileName,validNodeTypes,SpecialNodeType,uniformCSAreas);
validNodeTypes = 1:numel(tree);
treeVisualizer(tree,highlightedNodes,initiationSite{cellnumber},SpecialNodeType,newFigure,colorString,validNodeTypes,pixelUnits)

