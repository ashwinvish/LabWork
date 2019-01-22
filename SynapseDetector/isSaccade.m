function [ logic , type] =  isSaccade(cellID)
% Check is the cell in cellID is a saccadic cell
% logic = 1 if saccadic cell
% logic = 0 if not a saccadic cell
% type is the type of saccadic cell

sparseSaccadic = [76540 76622 76626 76697 76748 76750 76751 77122 77151 77238 ...
    77239 77240 77241 77437 77645 77708 77740 77826 78351 78541 78543 78545 ...
    78558 78572 78641 76611 77630 79058 79064 79067 79069 79080 79086 77374 ];

bushySaccadic = [76618 76625 76627 77132 77162 77163 77329 77434 77447 77460 ...
    77467 77797 77805 77848 78357 78358 79054 79055 79059 79062 79074 79077 ...
    79085 78544 78650 77390 77621 77636 77651 77656 ];

putativeBushySaccadic = [77342 77352 77336 77354 77373 78346 77433 77435 77453 77461 ];

lateralDSaccadic = [76629 76667 80185 76691 76692 77357 77389 77689 77806 77816 78601 77667 77684  ];

lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682 ];

unknownSaccadic = [78391 78434 78438 78440 78447 78448 78450 78540 78583 78649 79017 76878 77450 77670 77683 79045 ];

if ismember(cellID,sparseSaccadic)
    logic = 1;
    type = char('sparseSaccadic');
    
elseif ismember(cellID,bushySaccadic);
    logic = 1;
    type = char('bushySaccadic');
elseif ismember(cellID,putativeBushySaccadic)
    logic = 1;
    type = char('putativeBushySaccadic');
elseif ismember(cellID, lateralDSaccadic);
    logic= 1;
    type = char('lateralDSaccadic');
elseif ismember(cellID, lateralVSaccadic);
    logic = 1;
    type = char('lateralVSaccadic');
elseif ismember(cellID, unknownSaccadic)
    logic = 1;
    type = char('unknownSaccadic');
else
    logic = 0;
end

end