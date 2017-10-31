addpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts/jsonlab');
data = loadjson('77605.txt','SimplifyCell',1);
xStart = str2num(data.fragments(:,9:9+3))
xEnd = str2num(data.fragments(:,14:14+3))
yStart = str2num(data.fragments(:,19:19+3))
yEnd = str2num(data.fragments(:,24:24+3))
zStart = str2num(data.fragments(:,29:29+4))
zEnd = str2num(data.fragments(:,35:35+4))

xRange = [min(xStart) max(xEnd)]
yRange = [min(yStart) max(yEnd)]
zRange = [min(zStart) max(zEnd)]


