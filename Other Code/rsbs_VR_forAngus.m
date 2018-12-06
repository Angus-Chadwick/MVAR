%function rsbs_VR_forAngus

clear all

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

NID = 9;

%for NID=1:length(idinf),  % loop over sessions/blocks
    
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
    
    if 0,  % run analysis
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
        
        fprintf('%s  %2.0d:%s,B%d, N celllab=%d(of %d)\n',datestr(now),NID,idinf{NID}.id,idinf{NID}.block,length(find(~isnan(CELLLAB))),length(CELLLAB));
              
        % This part plots and then saves the calcium traces for each cell
        % of each type
        
        if 1, % plot traces
            for lab=1:length(TSK.(nam).TL.labB),  % loop over cell types
                
                % cell type info and indexing
                
                cnam = TSK.(nam).TL.labB{lab}; % cell type name
                if isempty(cnam), cnam='PYR'; end
                labix=find(CELLLAB==lab);  % find set of cell indices of this class 
                
                % set up figure properties
                
                interactive=1;toplot=0;              
                ext='-dpsc';fext='.ps';
                %ext ='-dtiff';fext    ='.tiff';
                filename = sprintf('%s%s_B%d_%s%s',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block,cnam,fext);               
                if toplot, delete(filename); end
                if interactive, vis='on'; else vis='off'; end
                figure('Visible',vis,'Units','centimeters','Position',[2,2,22,28]);
                six=1;subix=[5,3];
                
                % get data to plot
                
                for labix2=1:length(labix),  % go through cells in this class one by one
                                                           
                    layer = POOLCV(labix(labix2),1);  % find layer for this cell
                    B     = POOLCV(labix(labix2),2);  % further cell data
                    chn   = find(out.CRSPR{layer}.conversion(:,2)==B);  % find channel number for this cell
                    curtim=out.frmtim(out.sel{layer});  % extracts sampling times for the given layer
                    tix=find(~isnan(out.dFj{1}{layer}(chn,:)),1,'first'); % get time index of first non-NaN calcium signal
                                                           
                    subplot(subix(1),subix(2),six);six=six+1;
                    if ~isempty(tix),
                        strt=curtim(tix);stop=strt+1000;  % start time
                        trng=find(curtim>strt&curtim<stop); % time range
                        plot(curtim(trng),out.dFj{1}{layer}(chn,trng));
                        axis tight;
                    end
                    
                    % save plots
                    
                    title(sprintf('chn%d',chn),'FontSize',6,'Interpreter','none');
                    if six>prod(subix),
                        set(gcf,'PaperPositionMode','auto');
                        if interactive, pause; end
                        if toplot, eval(sprintf('print(''-f%d'',''-append'',''%s'',''%s'');',gcf,ext,filename)); end
                        clf;
                        six=1;
                    end
                end
                if six<=prod(subix),
                    set(gcf,'PaperPositionMode','auto');
                    if interactive, pause; end
                    if toplot, eval(sprintf('print(''-f%d'',''-append'',''%s'',''%s'');',gcf,ext,filename)); end
                end
                close all;
            end % for lab=1:TSK.(nam).TL.labB,
        end % if 0,
        
    end % if 1,
    
    
%--------------------------------------------------------------
% do some other analyses
    

    if 1, % run analysis
        clear ADDTRG;
        if 1, % rsbs_extract_data_trigger
            dbstop if error;
            ADDTRG=rsbs_extract_data_trigger(out,ADDINF,idinf,NID,EYE,EYEINF,TSK.(nam),RIXSES);          
            fname=sprintf('%sADDTRG2_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
            save(fname,'ADDTRG','-v7.3');
        end
        
        % population activity
        if 0, % rsbs_extract_data_trigger_model
            dbstop if error;
            ADDTRGpopact_model=rsbs_extract_data_trigger_popact_model(out,ADDINF,idinf,NID,EYE,EYEINF,TSK.(nam),RIXSES);
            
            if 1,
                CORLAG=ADDTRGpopact_model.CORLAG;
                fname=sprintf('%sADDTRGpopact_model_CORLAG_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
                save(fname,'CORLAG','-v7.3');
            end    
        end
        
        % single cell activity
        if 0, % rsbs_extract_data_trigger_modelcell_simple try fit model
            dbstop if error;
            ADDTRGcell_model=rsbs_extract_data_trigger_cell_model(out,ADDINF,idinf,NID,EYE,EYEINF,TSK.(nam),RIXSES,ADDTRG);
            
            CORLAG=ADDTRGcell_model.CORLAG;
            fname=sprintf('%sADDTRGcell_model_CORLAG_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
            save(fname,'CORLAG','-v7.3');
            %rsbs_predict_resp_corrlag
        end % if 0, % run analysis
        
    end % if 0, % run analysis
    
%end % for NID=1:length(idinf),
