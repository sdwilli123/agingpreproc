for rr=1:5
load(['/projectnb/fastfmri/sdwilli/aging/ag106a/behav/run' num2str(rr) '.mat']);

ti = allt-allt(1); 
i=2; ons = [] ; 
%get actual on times from behavioral files 
    while i < length(allt)
        if y(i-1)==1 && y(i) >1
            ons(end+1) = ti(i);
            i = i + (20*100); % skip to next chunk where an on might happen - at least 20 sec away x refresh rate of 100
        end 
        i=i+1;
    end 
%%
ons = ons;
save([ 'run' num2str(rr) '_onset.mat'], 'ons') 
end 