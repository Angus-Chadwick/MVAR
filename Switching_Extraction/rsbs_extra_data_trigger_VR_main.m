function [SIG]=rsbs_extra_data_trigger_VR_main(out,EYE, sel,BTRGLABEL,MOT)


%% This script takes input from rsbs_VR_forAngus.m, calculates data around specific trigger points on individual trials, and outputs this data as SIG. (AC)

%--------------------------------------------------------------------------
% edit rsbs_VR_main

frmtim = out.frmtim;

SFi = out.SFiLayer;  % sampling rate at a given layer (AC)


% cut out data around trigger
frmstrt = round(sel.start*SFi);  
frmstop = round((sel.lngth+sel.start)*SFi);

chnix=1;
cfginput.ero=0;

%% 

clear DAT DATTIM;

for layer=1:length(out.dFj{chnix}),
    
    DAT{layer}   = single(NaN([length(sel.trl),size(out.dFj{chnix}{layer},1),length([frmstrt:frmstop])]));  % matrix of trial x cell x time point around trigger (AC)
    DATTIM{layer}= single(NaN([length(sel.trl),length([frmstrt:frmstop])]));   % initialise window around trigger. sel.trl is probably time of trial? So trial x time point (AC)
   
    for T=1:length(sel.trl),
        
        [I]=find(frmtim(out.sel{layer})>sel.trl(T),1,'first');  % find index of first time point after trigger (AC)
        
        if isempty(I),I=length(frmtim(out.sel{layer})); end  % set to end value for ill-conditioned case (AC)
        
        rngix = [I+frmstrt:I+frmstop];  % indices for time vector around trigger for this trial (AC)
        
        if any((rngix<1)|(rngix>length(frmtim(out.sel{layer})))),  % if any of the range fall outwith the sampling vector for this layer (AC)
            
            ixok=(rngix>=1)&(rngix<=length(frmtim(out.sel{layer})));  % get subset which fall within the sampling vector (AC)
            rngix1=rngix(ixok); 
            rngix2=[1:length([frmstrt:frmstop])];rngix2=rngix2(ixok);
            
            if cfginput.ero,
                DAT{layer}(T,:,rngix2)= out.dFje{chnix}{layer}(:,rngix1);
            else
                DAT{layer}(T,:,rngix2)= out.dFj{chnix}{layer}(:,rngix1);  % get calcium data for this trial (AC)
            end
            DATTIM{layer}(T,rngix2)=frmtim(out.sel{layer}(rngix1));  % get time data for this trial (AC)
            
        else  % or if the full range is within the sampling vector
            
            if cfginput.ero,
                DAT{layer}(T,:,:)= out.dFje{chnix}{layer}(:,rngix);  
            else
                DAT{layer}(T,:,:)= out.dFj{chnix}{layer}(:,rngix);
            end
            DATTIM{layer}(T,:)=frmtim(out.sel{layer}(rngix)); % frmtim(rngix);
        end
    end
end % for layer=1:length(out.DAT{chnix}),

tax = [frmstrt:frmstop]./SFi; 

% compute PSTH's
clear SIG;
SIG.DAT=DAT;
SIG.MOT = MOT;
SIG.EYE = EYE;
SIG.sel=sel;
SIG.SFi=SFi;
SIG.CRSPR  = out.CRSPR;
SIG.BTRGLABEL = BTRGLABEL;
SIG.tax=tax;