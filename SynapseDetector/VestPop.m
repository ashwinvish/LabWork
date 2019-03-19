% Inputs and Outputs of Vestiblar Cells

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


 vestibularCellIds = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
     77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
     80672,80622,80572,80569,81122,78426];
 
 MVNs = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
     77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
     80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 ];
 
 allVest = [vestibularCellIds,MVNs];
 
 for i = 1:numel(vestibularCellIds)
     if isExistReRoot(vestibularCellIds(i)) == true
         VestCells(i) = InputsByClass(vestibularCellIds(i),df,1);
     end
 end
 
 
 for i = 1:numel(vestibularCellIds)
     % inputs
     numberOf.SaccadicInputs(i) = sum(VestCells(i).isSaccadic);
     numberOf.IntegratorInputs(i) = sum(VestCells(i).isIntegrator);
     numberOf.ContraInputs(i) = sum(VestCells(i).isContra);
     numberOf.RemainingInputs(i) = sum(VestCells(i).isEverythingElse);
     
      % output

      numberOf.SaccadicOutputs(i) = sum(isSaccade(VestCells(i).Outputs));
      numberOf.IntegratorOutputs(i) = sum(isIntegrator(VestCells(i).Outputs));
      numberOf.MotorOutputs(i) = sum(sum(isPostSynapseMotor(VestCells(i).Outputs)));
      numberOf.RemainingOutputs(i) = length(VestCells(i).Outputs)-(numberOf.SaccadicOutputs(i)+numberOf.IntegratorOutputs(i)+numberOf.MotorOutputs(i));
 end
 
 % plot location of the synapses
 SaccadicInputSynapses = vertcat(VestCells.PreSynCoordsTransformed);
 [a,b] = hist3([SaccadicInputSynapses(:,1),SaccadicInputSynapses(:,2)],[51 51]);
 imagesc(b{:},a');
 colormap(colorcet('L16','reverse',1));
 
 