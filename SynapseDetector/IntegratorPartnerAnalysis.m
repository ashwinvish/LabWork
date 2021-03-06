clc;
clear;

addpath(genpath('/Users/ashwin/Documents/'));

colors = cbrewer('qual','Dark2',10);
startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/11252018.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

confirmedALX = [76181 76201 76187 76184 76192 76197 ];
putativeALX = [76657 76623 76701 76661 76673 76690 76677 76624 77327 77580 ...
    77341 77344 77868 77872 77349 76634 77592 77578 77607 78629 77581 77591 ...
    77602 77821 77844 77845 77822 78406 78667 79022 79033 78404 78441 78421 79044 79046 79048 ];

ALXCellIDs = [confirmedALX,putativeALX];

confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
putativeDBX = [77582 77588 77595 77597 77598 77599 77605 79040 76399 76828 ...
    76829 78435 76826 76874 76289 76542 76832 76836 76838 76877 76883 76884 ...
    78400 78429 76383 76867 77594 ];

DBXCellIds  = [confirmedDBX,putativeDBX];

confirmedBARHL = [76198 76190 76193 76194 76195 76196 ];
putativeBARHL = [78452 80224];

BARHLCellIds = [confirmedBARHL,putativeBARHL];

%%

for i = 1:numel(ALXCellIDs)
    ALX(i) = InputsByClass(ALXCellIDs(i),df);
end

for i = 1:numel(DBXCellIds)
    DBX(i) = InputsByClass(DBXCellIds(i),df);
end

for i = 1:numel(BARHLCellIds)
    BARHL(i) = InputsByClass(BARHLCellIds(i),df);
end

%% plot the profile of input neurons

