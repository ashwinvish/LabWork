 clear

df = readtable('/Users/ashwin/Documents/SynapseDetector/11252018.csv');

triangle = [76749 76752 77374 77446 77455 78679 ];

sparseSaccadic = [76540 76622 76626 76697 76748 76750 76751 77122 77151 ...
    77238 77239 77240 77241 77437 77645 77708 77740 77826 78351 78545 78558 ...
    78572 78641 76611 77630 77374 80763 80850 80821 80801 80743 81007 80974 ...
    81002 80216 ];

bushySaccadicMedial = [76618 76625 76627 77132 77162 77163 77329 77434 77447 77460 ...
    77467 77797 77805 77848 78357 78358 79054 79059 78544 78650 77390 77621 ...
    77636 77651 77656 80995 ];

bushySaccadicLateral = [79067 79058 79069 79072 79055 79074 79062 79077 79085 ...
    79080 79086 79064 78576 78540 77146 78646 78542 78546 78541 78543 80728]; 

putativeBushySaccadic = [77342 77352 77336 77354 77373 78346 77433 77435 77453 77461];

lateralDSaccadic = [76629 76667 80185 76691 76692 77357 77389 77689 77806 ...
                    77816 78601 77667 77684 80542 ];

lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];

unknownSaccadic = [78583 78649 77670 80217 80746 80679 80804 80757 80647 ...
    80947 80939 81027 80943 80315 ];

allSaccade = [triangle,sparseSaccadic,bushySaccadicMedial,bushySaccadicLateral,...
    putativeBushySaccadic,lateralDSaccadic,lateralVSaccadic,unknownSaccadic];


%%

for i = 1:numel(allSaccade)
    motorDiff(i,:) = isMotor(allSaccade(i),df);
end
    
ABDheavy = allSaccade(find((motorDiff(:,2)+motorDiff(:,3)) > (motorDiff(:,4)+motorDiff(:,5))));
ABDiheavy = allSaccade(find((motorDiff(:,2)+motorDiff(:,3)) < (motorDiff(:,4)+motorDiff(:,5))));

%%



