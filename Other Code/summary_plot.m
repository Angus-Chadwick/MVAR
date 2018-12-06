%% This script extracts data and plots some basic information



% This section finds all sessions and blocks in the data folder and stores
% the file names in idinf

ROOT = '/nfs/nhome/live/angus/Documents/Interneuron_Data/';

RSG.savedirforAngus = [ROOT,'Saved_Data/'];

list=dir([RSG.savedirforAngus,'/*_dat*']);
clear idinf;
for NID=1:length(list),
    ix1=strfind(list(NID).name,'_B');  % find unique sessions, day and block number
    ix2=strfind(list(NID).name,'_dat');
    list(NID).name(1:ix1-1);
    idinf{NID}.id    = list(NID).name(1:ix1-1);
    idinf{NID}.block = str2num(list(NID).name(ix1+2:ix2-1));
end

%--------------------------------------------------------------------------

% This section extracts and plots different quantities of interest from the
% raw data files

clc;

if ~exist('TSK'), TSK=[]; end

if ~isfield(TSK,'SW'),
    TSK.SW=load([RSG.savedirforAngus,'/SW_TLTPDPOP.mat']);  % load in metadata for cells and animals
end

NID = 1;


    % get additional information about this session from the extinf file
    
    clear EXTINF;  
    fname=sprintf('%s%s_B%d_extinf',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
    load(fname,'EXTINF');
    nam    = EXTINF.nam;
    if strcmp(nam,'LN'), idinf{NID}.type='SD'; elseif strcmp(nam,'SW'), idinf{NID}.type='SWITCH'; end % determines if learning or switching
    RIXSES = EXTINF.RIXSES;  % used later for cell type classification
    EYE    = EXTINF.EYE;    
    EYEINF = EXTINF.EYEINF;
    
    clear out ADDINF;
    fname=sprintf('%s%s_B%d_dat',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
    load(fname,'out','ADDINF');
    
    
    
    %--------------------------------------------------------------
    % plot traces seperately for different interneuron type
    
    % This part seems to find the cell types and store them numerically as
    % CELLLAB. The numerical ordering corresponds to the entries in
    % TSK.SW.TL.labB
    

        POOLCV=[];
        for layer=1:4,
            dum=out.CRSPR{layer}.conversion(:,2);
            POOLCV  = cat(1,POOLCV,[repmat(layer,[size(dum,1),1]),dum]);
        end % for layer=1:length(TOT{ses}.SIG.DAT),
        
        % find cell type of each cell
        CELLLAB=NaN(size(POOLCV,1),1); % index into TL
        for chn=1:size(POOLCV,1),
            ix=find( (TSK.(nam).TL.rix==RIXSES(1))&(TSK.(nam).TL.ses==RIXSES(2))...
                &(TSK.(nam).TL.lay==POOLCV(chn,1))&(TSK.(nam).TL.B==POOLCV(chn,2)) );  
            if not(isempty(ix)),
                CELLLAB(chn) = TSK.(nam).TL.labJ(ix);
            else
                %warning(sprintf('cannot find chn %d',chn));
            end
        end
                    
        
        
        
        % This part plots and then saves the calcium traces for each cell
        % of each type
    
        
               
 for lab=1:length(TSK.(nam).TL.labB),  % loop over cell types
                  
     subplot(length(TSK.(nam).TL.labB), 1, lab);
     hold on
                
     % cell type info and indexing
                
     cnam = TSK.(nam).TL.labB{lab}; % cell type name
     if isempty(cnam), cnam='PYR'; end
     labix=find(CELLLAB==lab);  % find set of cell indices of this class 
                
     % get data to plot
                
     if lab == 1
         
          layer = POOLCV(labix,1);  % cell layers
          B     = POOLCV(labix,2);  % further cell data
          
          for LYR=1:4  
              
              layerix = find(layer==LYR);
              
              for CELL=1:length(layerix)
              
                  chn(CELL) = find(out.CRSPR{LYR}.conversion(:,2) == B(layerix(CELL)));  
                  
              end
              
              curtim=out.frmtim(out.sel{LYR});

              plot(curtim, mean(out.dFj{1}{LYR}(chn,:), 1))
              
              clear chn
              
          end
          
     else
              
     for labix2=1:length(labix),  % go through cells in this class one by one
                                                           
                    layer = POOLCV(labix(labix2),1);  % find layer for this cell
                    B     = POOLCV(labix(labix2),2);  % further cell data
                    chn   = find(out.CRSPR{layer}.conversion(:,2)==B);  % find channel number for this cell
                    curtim=out.frmtim(out.sel{layer});  % How does this extract the time? Is the time layer-dependent due to scanning properties?
                                                        
                    plot(curtim,out.dFj{1}{layer}(chn,:) + 2 * labix2);
                    axis tight;
                  
               
     end
     end
 end
    
        
