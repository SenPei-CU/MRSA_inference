function SILIinference()
% mex agentbasedsimulation.c
% mex agentbasedprob.c
load inward %in which ward patient stay
load wardcap %ward capacity
tmstep=1;%observation at every tmstep week
Nmax=743599;%total number of patients
%parameters: alpha and p used in the synthetic modelling
alphal=1.0/365; alphau=1.0/175;%decolonization rate
pkl=0.15; pku=0.25;%observation rate = pxalpha
load observation %individual-level observation: patient id, diagnose week
load paratruth %parameters used in generating synthetic outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_times=52;%test model in 52 weeks
Cl=0; Cu=0.05;%initial colonization rate for patients in hospital
num_ens=100;%number of ensemble
appear=zeros(Nmax,1);%record if nodes appear before current time
x=xtruth(1:5);%parameters used in generating synthetic outbreak
%beta0,I0=0 (imporation rate of infection),
%C0=gamma (importation rate of colonization),theta,D
x=x*ones(1,num_ens);
para_ens=zeros(Nmax,3);
para_ens(:,1)=(alphau+alphal)/2;%decolonization rate
para_ens(:,2)=(pku+pkl)/2;%observation rate
para_ens(:,3)=para_ens(:,1);%decolonization rate of observed colonization
%initialize
Sl=1-Cu; Su=1-Cl;
%load the first network
burnin=52;%synthetic outbreak starts from week 53
filePattern = fullfile(pwd,'net/net');
load([filePattern,num2str(burnin+1),'.mat']);%start from week 25
existingnodes=unique(nl(nl(:,2)>0,1));
appear(existingnodes)=1;
S_cnt=zeros(length(existingnodes),num_ens+1);
C_cnt=zeros(length(existingnodes),num_ens+1);

S_cnt(:,1)=existingnodes;
S_cnt(:,2:num_ens+1)=Sl+rand(length(existingnodes),num_ens)*(Su-Sl);
C_cnt(:,1)=existingnodes;
C_cnt(:,2:num_ens+1)=1-S_cnt(:,2:num_ens+1);

wardnum=1041;
wardcontam=zeros(wardnum,num_ens);%contamination rate in each ward
wardcontam_post=zeros(wardnum,num_ens,num_times);
infnodes=diagnoserec(:,1);%infected nodes
%%%%%%%%%%%%%%posterior colonization probability for observed colonization
Cpost_rec=NaN(size(diagnoserec,1),num_ens+1,num_times);
appear_last=zeros(Nmax,1);%nodes appeared until last time step
%%%%%%%%%%%%%%%start looping at each day
datetime
for t=1:num_times
    t
    %%%%%%%%%%%%%calculate Cpost_obs, back propagation, at the beginning
    %%%%%%%%%%%%%of week t
    %Cprior_obs: prior colonization probability for observed carriers
    %Cpost_obs: posterior colonization probability for observed carriers
    newnodes=setdiff(find(appear>0),find(appear_last>0));%find newly admitted nodes
    [Cpost_obs,Cprior_obs]=getCpost_obs(x,S_cnt,C_cnt,appear,t,para_ens,infnodes,wardcontam,diagnoserec,burnin,inward,wardcap,newnodes);
    for i=1:size(Cpost_obs,1)
        inode=Cpost_obs(i,1);
        tempid=find(diagnoserec(:,1)==inode);
        for j=1:length(tempid)
            Cpost_rec(tempid(j),:,t)=Cpost_obs(i,:);
        end
    end
    %%%%%%%%%%%%%%update inode and neighbors
    cntinfnum=size(Cpost_obs,1);%number of infected nodes in the network
    for i=1:cntinfnum
        inode=Cpost_obs(i,1);
        inodeindex = find(C_cnt(:,1)==inode);
        %Cpost_obs dim: cntinfnum,num_ens+1
        if (sum(isnan(Cpost_obs(i,2:end)))==0)&&(sum(isnan(C_cnt(inodeindex,2:end)))==0)
            %%%%%%%get prior
            obs_ens = C_cnt(inodeindex,2:end);
            prior_var = var(obs_ens);
            if (prior_var==0)||(isnan(prior_var))
                obs_ens=obs_ens.*(1+0.01*randn(1,num_ens));
                obs_ens(obs_ens>=1)=0.999+0.001*rand(1,sum(obs_ens>=1));
                obs_ens(obs_ens<=0)=0.001*rand(1,sum(obs_ens<=0));
                prior_var=var(obs_ens);
            end
            %%%%update the infected nodes
            %C
            temp=Cpost_obs(i,2:end);
            %avoid 0 or 1
            temp(temp>=1)=0.999+0.001*rand(1,sum(temp>=1));
            temp(temp<=0)=0.001*rand(1,sum(temp<=0));
            dy=temp-C_cnt(inodeindex,2:end);
            C_cnt(inodeindex,2:end)=temp;
            %S
            S_cnt(inodeindex,2:end)=1-C_cnt(inodeindex,2:end);
            %%%%loop through all neighbors
            for j=part(inode):part(inode+1)-1
                nei=nl(j,1);
                neiindex=find(C_cnt(:,1)==nei);
                if ~isempty(neiindex)%neighbor exists at the beginning of week t
                    %C
                    r=cov(obs_ens,C_cnt(neiindex,2:end));
                    r=r(2,1)/prior_var;
                    dx=r*dy;
                    temp=C_cnt(neiindex,2:end)+dx;
                    %avoid 0 or 1
                    temp(temp>=1)=0.999+0.001*rand(1,sum(temp>=1));
                    temp(temp<=0)=0.001*rand(1,sum(temp<=0));
                    C_cnt(neiindex,2:end)=temp;
                    %S
                    S_cnt(neiindex,2:end)=1-C_cnt(neiindex,2:end);
                end
            end
        end
    end
    %%%%%%%%%%%%%%%step forward for one step
    filePattern = fullfile(pwd,'net/net');
    load([filePattern,num2str(burnin+t),'.mat']);
    appear_last=appear;
    [S_cnt,C_cnt,appear,wardcontam]=modelprob(x,S_cnt,C_cnt,appear,burnin+t,tmstep,para_ens,nl,part,deg,wardcontam,inward,wardcap);
    wardcontam_post(:,:,t)=wardcontam;
    
    patients=unique(nl(nl(:,8)>0,1));%patients in hospital on the last day of the week
    C_infer=C_cnt(ismember(C_cnt(:,1),patients),:);%inferred colonization probability
end
datetime
save SILIinference.mat C_infer

%plot ROC curve
load state_truth
load('SILIinference.mat', 'C_infer')
truth=state_truth(:,end);%true states
colonized=find(truth==1);
patients=C_infer(:,1);
inhospcolonized=zeros(Nmax,1);
inhospcolonized(intersect(colonized,patients))=1;
rank=zeros(Nmax,2);
rank(:,1)=(1:Nmax)';
rank(C_infer(:,1),2)=mean(C_infer(:,2:end),2);
rank=rank(patients,:);
rank=sortrows(rank,-2);
labels=inhospcolonized(rank(:,1));
scores=rank(:,2);
[X,Y,~,~] = perfcurve(labels,scores,1);
plot(X,Y,'LineWidth',2);
xlabel('False positive rate')
ylabel('True positive rate')

function [Cpost_obs,Cprior_obs]=getCpost_obs(x,S_cnt,C_cnt,appear,Tstart,para_ens,infnodes,wardcontam,diagnoserec,burnin,inward,wardcap,newnodes)
num_ens=size(S_cnt,2)-1;
tmstep=1;
%S_cnt and C_cnt give the prior
%run through all infected nodes in network to get likelihood
inodes=intersect(infnodes,S_cnt(:,1));%find infected nodes in network
inodesafter=diagnoserec(diagnoserec(:,2)>=Tstart,1);%find infected nodes infected after Tstart
inodes=intersect(inodes,inodesafter);
inodes=intersect(inodes,newnodes);%infected nodes that firs appear
num_inodes=length(inodes);
Cpost_obs=zeros(num_inodes,num_ens+1);
Cpost_obs(:,1)=inodes;
Cprior_obs=Cpost_obs;
% [Tstart,num_inodes]
tic
%%%%%%%%%prepare ensemble for all infected nodes
xtemp=mean(x,2)*ones(1,2*num_inodes);
%C for all infected nodes + S for all infected nodes
S_cnt_temp=zeros(size(S_cnt,1),1+2*num_inodes); C_cnt_temp=zeros(size(C_cnt,1),1+2*num_inodes);
S_cnt_temp(:,1)=S_cnt(:,1); C_cnt_temp(:,1)=C_cnt(:,1);
%%%calculate mean probability
S_cnt_temp(:,2:end)=mean(S_cnt(:,2:end),2)*ones(1,2*num_inodes);
C_cnt_temp(:,2:end)=mean(C_cnt(:,2:end),2)*ones(1,2*num_inodes);
for i=1:num_inodes
    inode=inodes(i);
    index=find(S_cnt(:,1)==inode);
    c_col=1+i;
    S_cnt_temp(index,c_col)=0; C_cnt_temp(index,c_col)=1;%set P(Xt=C)=1
    s_col=1+num_inodes+i;
    S_cnt_temp(index,s_col)=1; C_cnt_temp(index,s_col)=0;%set P(Xt=S)=1
end
appear1=appear;
wardcontam1=mean(wardcontam,2)*ones(1,2*num_inodes);
%find later infection times
itimes=unique(diagnoserec(diagnoserec(:,2)>=Tstart,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P(Obs|Xt=C)
lh_c=ones(1,num_inodes);%likelihood
%P(Obs|Xt=S)
lh_s=ones(1,num_inodes);%likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%run until each itimes
for ti=1:length(itimes)
%     [ti,length(itimes),num_inodes]
    %find the start and end time of each interval
    if ti==1
        tstart=Tstart;
    else
        tstart=itimes(ti-1)+1;
    end
    tend=itimes(ti);
    for t=tstart:tend
        filePattern = fullfile(pwd,'net/net');
        load([filePattern,num2str(burnin+t),'.mat']);
        [S_cnt_temp,C_cnt_temp,appear1,wardcontam1]=modelprob(xtemp,S_cnt_temp,C_cnt_temp,appear1,burnin+t,tmstep,para_ens,nl,part,deg,wardcontam1,inward,wardcap);
    end
    %%%%%calculate likelihood
    %%%find infected nodes at itimes(ti)
    cntinodes=diagnoserec(diagnoserec(:,2)==itimes(ti),1);
    for j=1:length(cntinodes)
        cntinode=cntinodes(j);
        %%%%%%%%%%%%for condition P(Xt=C)=1
        Pc=C_cnt_temp(C_cnt_temp(:,1)==cntinode,2:1+num_inodes);
        %calculate likelihood, chain
        if (~isempty(Pc))&&(sum(isnan(Pc))==0)
            lh_c=lh_c.*Pc;
        end
        %%%%%%%%%%%%for condition P(Xt=S)=1
        Pc=C_cnt_temp(C_cnt_temp(:,1)==cntinode,1+num_inodes+1:end);
        %calculate likelihood, chain
        if (~isempty(Pc))&&(sum(isnan(Pc))==0)
            lh_s=lh_s.*Pc;
        end
        %%%%%%%%%%set P(cntinode=C)=1
        C_cnt_temp(C_cnt_temp(:,1)==cntinode,2:end)=1;
        S_cnt_temp(S_cnt_temp(:,1)==cntinode,2:end)=0;
    end
end
for i=1:num_inodes
    inode=inodes(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %normalize
    index=find(S_cnt(:,1)==inode);
    prior_s=S_cnt(index,2:end); prior_c=C_cnt(index,2:end);
    Cprior_obs(i,2:end)=prior_c;
    Cpost_obs(i,2:end)=(prior_c*lh_c(i))./(prior_c*lh_c(i)+prior_s*lh_s(i)); 
end
toc
