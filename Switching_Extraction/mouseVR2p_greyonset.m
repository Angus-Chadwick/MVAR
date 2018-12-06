function STMONSET = mouseVR2p_runonset(out,AUXPOS,idinf,PARAM);

if ismember(idinf.type,{'SD','SWITCH'}),
    if strcmp(idinf.type,'SD'),
        %------------------------------------------------------------------
        % make triggers for circle onset and grey onset
        % edit rsbs_VR_speed.m
        clear onsetcircle onsetgrey;
        for cnd=1:2,
            if cnd==1, dum=PARAM.gratingonnormal; end
            if cnd==2, dum=PARAM.gratingon2; end
            onsetcircle{cnd} = NaN(size(dum));
            onsetgrey{cnd}   = NaN(size(dum));
            for i=1:length(dum), %i=2
                [Y,ix]=min(abs(AUXPOS.tax-dum(i)./1000));
                % AUXPOS.y(ix-1:ix+1)
                if AUXPOS.y(ix-1)~=0&AUXPOS.y(ix)==0,
                    error('here');
                end
                % find onset of preceding circle corridor
                ix2=find(AUXPOS.y(1:ix-1)~=0,1,'last')+1;
                
                % find onset of grey corridor
                ix3=ix2+find(AUXPOS.x(ix2:ix-1)>0,1,'first')-1;
                
                if 0,
                    AUXPOS.y(ix-2:ix+2)'
                    AUXPOS.y(ix2-2:ix2+2)'
                    [AUXPOS.x(ix3-2:ix3+2)';AUXPOS.y(ix3-2:ix3+2)']
                end
                
                if ~isempty(ix2),onsetcircle{cnd}(i)=AUXPOS.tax(ix2).*1000;end
                if ~isempty(ix3),onsetgrey{cnd}(i)  =AUXPOS.tax(ix3).*1000;end
            end % for i=1:length(dum),
        end % for cnd=1:2,
        
        STMONSET.onsetcircle = onsetcircle;
        STMONSET.onsetgrey   = onsetgrey;
    end % if strcmp(idinf{NID},'SD'),
    
    if strcmp(idinf.type,'SWITCH'),
        runup=AUXPOS.y==0;
        runupix=find(diff(runup)==1)+1;
        %AUXPOS.y(runupix(1)-1:runupix(1)+1);
        if 0,
            figure;plot(AUXPOS.y);
            hold on;plot([runupix,runupix],get(gca,'YLim'),'k');
        end
        
        runupix=AUXPOS.tax(runupix);
        STMONSET.onsetgrey = runupix.*1000;
    end % if strcmp(idinf{NID},'SWITCH'),
else
    STMONSET=[];
end
