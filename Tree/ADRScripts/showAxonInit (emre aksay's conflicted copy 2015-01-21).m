function [ output_args ] = showAxonInit(fileName,cellnumber)
%showAxonInit -- highlight axon initiation site

% hardcoded initiation sites
initiationSite{4} = [39 46 47 48 52 49 50 51 53 54 55 56 57 58 59 60 62 63 64 65 66 61]; 
initiationSite{5} = [56 75 91 93 104 105 106 107 108 109 110 111 112 114 115 116 113 ];
initiationSite{6} = [11 17 20 26 82 95 96 97 98 133 138 143 147 149 150 156 154 148 151 152 122 125 127 132 128 130 134 140 142 153 165 139 141 155 157 159 160 163 162 158 164 166 131 117 119 120 135 145 161 167 168 169 170 171 172 129 136 144 146 137 56 66 68 72 75 83 117 86 99];
initiationSite{7} = [9 30 70 71 76 92 105 109 113 103 104 106 107 110 111 112];
initiationSite{21} = [1 91 34 52 80 82 83 86 88 90 103 104 107 111 112 113 114 115 92 81 84 87 89 95 96 98 99 100 93 94 101 102]; 
SpecialNodeType = []; 
validNodeTypes = [-1:5];
uniformCSAreas = [];

highlightedNodes = 1;
newFigure = true;
colorString = {'g'};
pixelUnits = false;

tree = generateIrreducibleDoubleLinkedTree(fileName{cellnumber},validNodeTypes,SpecialNodeType,uniformCSAreas);
treeVisualizer(tree,highlightedNodes,initiationSite{cellnumber},SpecialNodeType,newFigure,colorString,validNodeTypes,pixelUnits)


end

