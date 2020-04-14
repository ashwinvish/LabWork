% plot K values
colorSchemes;

fileID = fopen('/Users/ashwin/Documents/LabWork/SynapseDetector/slope_data/abd_inter_fraction.txt');
abd_inter = textscan(fileID,'%f');
fclose(fileID);

fileID = fopen('/Users/ashwin/Documents/LabWork/SynapseDetector/slope_data/bin_edges.txt');
bin_edges = textscan(fileID,'%f');
fclose(fileID);

bin_edges = -19.5:1:29.5;


fileID = fopen('/Users/ashwin/Documents/LabWork/SynapseDetector/slope_data/abd_motor_fraction.txt');
abd_motor = textscan(fileID,'%f');
fclose(fileID);

fileID = fopen('/Users/ashwin/Documents/LabWork/SynapseDetector/slope_data/integrator_fraction.txt');
integrator = textscan(fileID,'%f');
fclose(fileID);

fileID = fopen('/Users/ashwin/Documents/LabWork/SynapseDetector/slope_data/vestibular_fraction.txt');
vestibular = textscan(fileID,'%f');
fclose(fileID);

fileID = fopen('/Users/ashwin/Documents/LabWork/SynapseDetector/slope_data/average_ks.txt');
average_ks = textscan(fileID,'%f');
fclose(fileID);

figure;
subplot(4,4,1)
plot(-19.5:1:29.5,abd_motor{1},'color',ABDcolor,'LineWidth',2);
hold on;
plot(-19.5:1:29.5,abd_inter{1},'color',ABDicolor,'LineWidth',2);
plot(-19.5:1:29.5,integrator{1},'color',SaccABDcolor,'LineWidth',2);
plot(-19.5:1:29.5,vestibular{1},'color',[255,127,0]./255,'LineWidth',2);
box off;
axis square;
offsetAxes(gca);



