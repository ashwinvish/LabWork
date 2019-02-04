function [logic] = isVestibular(cellID, df)


 vestibularCellIds = [76631 76700 76656 76655 76660 76679 77393 77395 77375 ...
     77815 77803 77807 80591 80579 77255 77697 77364 80223 81126 80760 81117 ...
     80672 80622 80572 ];
 
logic = ismember(cellID, vestibularCellIds');

end
