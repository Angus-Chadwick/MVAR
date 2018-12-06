function RUNONSET = mouseVR2p_runonset(out,ADDINF,idinf,RUNONSET);


AUXPOS = ADDINF.AUXPOS;
%--------------------------------------------------------------------------
% detect running on and offset

S=AUXPOS;
S.xraw=medfilt1(S.xraw,30)*774;


r=ADDINF.valveonind'./1000;
ons             = RUNONSET.offs./1000; % offset of still = onset of running
offs            = RUNONSET.ons./1000;  % onset of still  = offset of running

if 0,
    % only on and offsets in runup corridor, and stay at least 1s in runup corridor
    dum=[RUNONSET.offsrunupexit-RUNONSET.offs]/1000;
    ons             = RUNONSET.offs(dum>1); % offset of still = onset of running
    offs            = RUNONSET.ons(dum>1);  % onset of still  = offset of running
    
       

  
end

if ismember(idinf.type,{'SD','SWITCH'}),
    %if strcmp(idinf.type,'SD'),
               
        onsok = zeros(size(ons));
        for i=1:length(r),
            [Y,ix1]=min(abs(AUXPOS.tax-r(i))); % reward
            
            % find when mouse exists corridor
            ix2=ix1+find(AUXPOS.y(ix1:end)==0,1,'first')-1;
            if isempty(ix2), ix2=length(AUXPOS.tax); end
            %diff(AUXPOS.tax([ix1,ix2]))
            
            % find any offsets that occur in between
            ix3=find(ons>AUXPOS.tax(ix1) & ons<AUXPOS.tax(ix2) ); %ix3=ix3(1);
            if ~isempty(ix3),
                %[i;ix3]
                onsok(ix3(1))=1;
            end
        end

    %end % if strcmp(idinf{NID},'SD'),
end
RUNONSET.ons           = ons.*1000;
RUNONSET.onsok         = onsok;
RUNONSET.onsRwd        = RUNONSET.ons( onsok==1 );


if 0,
    figure;plot(S.tax,Speed);
    hold on;plot(get(gca,'XLim'),[6,6],'k');
    hold on;plot([ons,ons],get(gca,'YLim'),'r');
    hold on;plot([offs,offs],get(gca,'YLim'),'g');

    g=ADDINF.gratingonnormal./1000;
    close all;
    figure;
    %plot(S.tax,S.xraw./max(S.xraw));
    plot(S.tax,S.xraw);
    %hold on;plot(S.tax,-S.y./max(-S.y),'k');
    hold on;plot([RUNONSET.onsRwd,RUNONSET.onsRwd]./1000,get(gca,'YLim'),'b','LineWidth',2);
    %hold on;plot([offs,offs],get(gca,'YLim'),'g','LineWidth',2);
    hold on;plot([r,r],get(gca,'YLim'),'r','LineWidth',2);
    %hold on;plot([g,g],get(gca,'YLim'),'b');
end

