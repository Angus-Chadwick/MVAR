function [SIG,SIGPOP]=rsbs_extra_data_trigger_behavdata_model(out,EYE, sel,BTRGLABEL,MOT,LN,RIXSES)

% see rsbs_popactivity
POOLCV=[];
for layer=1:4,
    dum=out.CRSPR{layer}.conversion(:,2);
    POOLCV  = cat(1,POOLCV,[repmat(layer,[size(dum,1),1]),dum]);
end % for layer=1:length(TOT{ses}.SIG.DAT),

out.POOLDAT=[];
for layer=1:4,
    DUM=interp1(out.frmtim(out.sel{layer}),out.dFj{1}{layer}',out.frmtim);
    out.POOLDAT=cat(1,out.POOLDAT,DUM');
end

SFi = out.SFi;

layer = 1;
out=rmfield(out,'sel');
out.sel{layer} = [1:length(out.frmtim)]';

% cut out data around trigger
frmstrt = round(sel.start*out.SFi);
frmstop = round((sel.lngth+sel.start)*out.SFi);

DUM = out.POOLDAT;

clear DAT DATTIM;
DAT{layer}   =single(NaN([length(sel.trl),size(DUM,1),length([frmstrt:frmstop])]));
DATTIM{layer}=single(NaN([length(sel.trl),length([frmstrt:frmstop])]));
for T=1:length(sel.trl),
    [I]=find(out.frmtim>sel.trl(T),1,'first');
    
    if isempty(I),I=length(out.frmtim); end
    rngix = [I+frmstrt:I+frmstop];
    if any((rngix<1)|(rngix>length(out.frmtim))),
        ixok=(rngix>=1)&(rngix<=length(out.frmtim));
        rngix1=rngix(ixok);
        rngix2=[1:length([frmstrt:frmstop])];rngix2=rngix2(ixok);
        DAT{layer}(T,:,rngix2)= DUM(:,rngix1);
        DATTIM{layer}(T,rngix2)=out.frmtim(rngix1);
    else
        DAT{layer}(T,:,:)= DUM(:,rngix);
        DATTIM{layer}(T,:)=out.frmtim(rngix); % out.frmtim(rngix);
    end
end
tax = [frmstrt:frmstop]./out.SFi;

% compute PSTH's
clear SIG;
SIG.DAT=DAT;
SIG.MOT = MOT;
SIG.EYE = EYE;
SIG.sel=sel;
SIG.SFi=SFi;
SIG.CRSPR  = out.CRSPR;
SIG.BTRGLABEL = BTRGLABEL;
SIG.POOLCV    = POOLCV;
SIG.tax=tax;


%--------------------------------------------------------------------------
% rsbs_popactivity

okchn = find(~isnan(nanmean(out.POOLDAT,2)));

CLASSLAB={'','PV','SOM','VIP'};

% find cell type of each cell
CELLLAB=NaN(size(out.POOLDAT,1),1); % index into TL
for chn=1:size(out.POOLDAT,1),
    ix=find( (LN.TL.rix==RIXSES(1))&(LN.TL.ses==RIXSES(2))...
        &(LN.TL.lay==POOLCV(chn,1))&(LN.TL.B==POOLCV(chn,2)) );;
    if not(isempty(ix)),
        CELLLAB(chn) = LN.TL.labJ(ix);
    else
        warning(sprintf('cannot find chn %d',chn));        
    end
end

CLASSNUM=zeros(size(CELLLAB));
for CLASS=1:length(CLASSLAB),
    CLASSNUM(CELLLAB==strmatch(CLASSLAB{CLASS},LN.TL.labB,'exact'))=CLASS;
end

% seperately for classes rsbs_medianavg_split
pc=NaN([length(CLASSLAB),size(out.POOLDAT,2)]); % pca
sm=NaN([length(CLASSLAB),size(out.POOLDAT,2)]); % summed activity
for C=1:length(CLASSLAB),
    ix=find(CLASSNUM==C);
    
    if ~isempty(ix),
        dat=out.POOLDAT(ix,:);
        okchnclass=find(~isnan(nanmean(dat,2))); % ok channels of this class
        oksmp=find(all(~isnan(dat(okchnclass,:)),1));
        
        [COEFF,SCORE] = princomp(dat(okchnclass,oksmp)');
        
        pc(C,oksmp)=SCORE(:,1)';
        sm(C,oksmp)=nansum(dat(okchnclass,oksmp),1);       
    end % if ~isempty(ix),
end
if 0,
    for C=1:length(CLASSLAB),
        close all;
        figure;%plot(out.frmtim,dat);
        subplot(2,1,1);hold on;plot(out.frmtim,sm(C,:),'Color',[1 0 0],'LineWidth',2);
        subplot(2,1,2);hold on;plot(out.frmtim,pc(C,:),'Color',[1 0 0],'LineWidth',2);
        
        pause;close;
    end
end

pc=(pc-repmat(nanmean(pc,2),[1,size(pc,2)]))./repmat(nanstd(pc,[],2),[1,size(pc,2)]);
sm=(sm-repmat(nanmean(sm,2),[1,size(sm,2)]))./repmat(nanstd(sm,[],2),[1,size(sm,2)]);
if 0,
    for C=1:length(CLASSLAB),
        close all;
        figure;%plot(out.frmtim,dat);
        hold on;plot(out.frmtim,sm(C,:),'Color',[1 0 0],'LineWidth',2);
        hold on;plot(out.frmtim,pc(C,:),'Color',[0 0 1],'LineWidth',2);
        
        pause;close;
    end
end
TRACE.pc = pc;
TRACE.sm = sm;

%--------------------------------------------------------------------------
% rsbs_extra_data_trigger_popact_VR_main

DUM = cat(1,TRACE.pc,TRACE.sm);

clear DAT DATTIM;
DAT{layer}   =single(NaN([length(sel.trl),size(DUM,1),length([frmstrt:frmstop])]));
DATTIM{layer}=single(NaN([length(sel.trl),length([frmstrt:frmstop])]));
for T=1:length(sel.trl),
    [I]=find(out.frmtim>sel.trl(T),1,'first');
    
    if isempty(I),I=length(out.frmtim); end
    rngix = [I+frmstrt:I+frmstop];
    if any((rngix<1)|(rngix>length(out.frmtim))),
        ixok=(rngix>=1)&(rngix<=length(out.frmtim));
        rngix1=rngix(ixok);
        rngix2=[1:length([frmstrt:frmstop])];rngix2=rngix2(ixok);
        DAT{layer}(T,:,rngix2)= DUM(:,rngix1);        
        DATTIM{layer}(T,rngix2)=out.frmtim(rngix1);
    else
        DAT{layer}(T,:,:)= DUM(:,rngix);
        DATTIM{layer}(T,:)=out.frmtim(rngix); % out.frmtim(rngix);
    end
end

% compute PSTH's
clear SIGPOP;
SIGPOP.DAT = DAT;
SIGPOP.MOT = MOT;
SIGPOP.EYE = EYE;
SIGPOP.sel = sel;
SIGPOP.SFi = out.SFi;
SIGPOP.CRSPR  = out.CRSPR;
SIGPOP.BTRGLABEL = BTRGLABEL;
SIGPOP.tax=tax;
SIGPOP.label = {'pcPYR','pcPVB','pcSOM','pcVIP','smPYR','smPVB','smSOM','smVIP'};