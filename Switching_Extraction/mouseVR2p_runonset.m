function RUNONSET = mouseVR2p_runonset(out,AUXPOS,idinf);

%--------------------------------------------------------------------------
% detect running on and offset

S=AUXPOS;
S.xraw=medfilt1(S.xraw,30);

% edit agk_rsbs_ROI_CRSP_rsp_plot_runningonsets.m

Speed = S.xraw*774;
still = Speed<6&Speed>-6;
ons=S.tax(find(diff(still)==1)+1);   % start of still
offs=S.tax(find(diff(still)==-1)+1); % end of still

if offs(1)<ons(1);offs(1) = [];end
if ons(end)>offs(end);ons(end) = [];end
sdur = offs-ons;
sdurrng = [.5  Inf];

trgsok = sdur>sdurrng(1)&sdur<=sdurrng(2);
trgsok(ons<out.frmtim(1)&offs>out.frmtim(end))=0;

ons=ons(trgsok);
offs=offs(trgsok);

if ismember(idinf.type,{'SD','SWITCH'}),
    if strcmp(idinf.type,'SD'),
        % restrict to onset and offset that both occur in circle runup corridor
        % and where mouse stays in runup corridor at least 0.5 s after
        % running onset
        offsrunupexit=NaN(size(offs));
        offsthreshold=NaN(size(offs));
        for i=1:length(ons),
            % for each stationary offset, find how long before exits runup
            % also find how long above stationary threshold
            [Y,ix1]=min(abs(AUXPOS.tax-ons(i)));
            [Y,ix2]=min(abs(AUXPOS.tax-offs(i)));
            if ((AUXPOS.y(ix1)==0)&(AUXPOS.x(ix1)>0))&...
                    ((AUXPOS.y(ix2)==0)&(AUXPOS.x(ix2)>0)),
                ix3=ix2+find(AUXPOS.y(ix2:end)~=0,1,'first')-2; % for each stationary offset, find how long before exits runup
                if isempty(ix3),ix3=length(AUXPOS.y);end
                ix4=ix2+find(still(ix2:end)==1,1,'first')-1;    % find how long above stationary threshold
                if isempty(ix4),ix4=length(AUXPOS.y);end
                offsrunupexit(i)=S.tax(ix3);
                offsthreshold(i)=S.tax(ix4);
                if 0,
                    AUXPOS.y(ix2:ix3)'
                    Speed(ix2:ix4)'
                end
            end
        end
    end % if strcmp(idinf{NID},'SD'),
    
    if strcmp(idinf.type,'SWITCH'),
        % restrict to onset and offset that both occur in runup corridor
        % and where mouse stays in runup corridor at least 0.5 s after
        % running onset
        offsrunupexit=NaN(size(offs));
        offsthreshold=NaN(size(offs));
        for i=1:length(ons), % i=14
            % for each stationary offset, find how long before exits runup
            % also find how long above stationary threshold
            [Y,ix1]=min(abs(AUXPOS.tax-ons(i)));
            [Y,ix2]=min(abs(AUXPOS.tax-offs(i)));
            
            % for switching
            if ((AUXPOS.y(ix1)==0)&(AUXPOS.x(ix1)>0)),
                ix3=ix2+find(AUXPOS.y(ix2:end)~=0,1,'first')-2; % for each stationary offset, find how long before exits runup
                if isempty(ix3),ix3=length(AUXPOS.y);end
                ix4=ix2+find(still(ix2:end)==1,1,'first')-1;    % find how long above stationary threshold
                if isempty(ix4),ix4=length(AUXPOS.y);end
                
                offsrunupexit(i)=S.tax(ix3);
                offsthreshold(i)=S.tax(ix4);
                if 0,
                    AUXPOS.y(ix2:ix3)'
                    Speed(ix2:ix4)'
                end
            end
        end
    end % if strcmp(idinf{NID},'SWITCH'),
    
    RUNONSET.offsrunupexit = offsrunupexit.*1000;
    RUNONSET.offsthreshold = offsthreshold.*1000;
end
RUNONSET.ons           = ons.*1000;
RUNONSET.offs          = offs.*1000;


if 0,
    figure;plot(S.tax,Speed);
    hold on;plot(get(gca,'XLim'),[6,6],'k');
    hold on;plot([ons,ons],get(gca,'YLim'),'r');
    hold on;plot([offs,offs],get(gca,'YLim'),'g');
end

