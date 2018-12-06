 %% This function gets cell class information and extracts calcium signals and time vectors for each cell


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