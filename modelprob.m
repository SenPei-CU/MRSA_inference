function [S_cnt,C_cnt,appear,wardcontam]=modelprob(x,S_cnt,C_cnt,appear,t,tmstep,para_ens,nl,part,deg,wardcontam,inward,wardcap)
%Run master equations for one week
%x:beta - transmission rate; I0 - new introductory infection; C0=gamma - new
%introductory colonization; theta - spillover rate; D - environmental clearance period
%%%%%%input
%S_ens - individual susceptible probability
%I_ens - individual colonized probability
%appaer - nodes ever appeared at t-1
%t - start week
%tmstep - simulation length (week)
%para_ens - disease parameters for each individual drawn from distributions
%nl, part, deg - contact network structure
%wardcontam - the environmental force of infection
%inward - in which ward patients stay
%wardcap - ward capacity
%%%%%%output
%S_ens - individual susceptible probability
%I_ens - individual colonized probability
%appear - nodes ever appeared at t
%wardcontam - updated environmental force of infection

num_ens=size(x,2);
%transform D to lambda
x(end,:)=ones(1,num_ens)./x(end,:);
Nmax=743599;
for t1=t:t+tmstep-1
    nodes=find(deg>0);
    wardcnt=inward{t1};
    oldnodes=find(appear>0);
    appear(deg>0)=1;
    appearnodes=find(appear>0);
    newnodes=setdiff(appearnodes,oldnodes);
    tempS_cnt=zeros(length(appearnodes),num_ens+1);
    tempC_cnt=zeros(length(appearnodes),num_ens+1);
    tempS_cnt(:,1)=appearnodes;
    tempC_cnt(:,1)=appearnodes;
    for k=1:num_ens
        C0=x(3,k);%importation rate of colonized patients
        para=x(1:5,k);%beta0,I0,C0,theta,lambda
        wardcontamcnt=wardcontam(:,k);
        alpha=para_ens(:,1);%decolonization rate
        pk=para_ens(:,2);%observation rate
        %the current state
        Scnt=zeros(Nmax,1);
        Scnt(S_cnt(:,1))=S_cnt(:,k+1);
        Ccnt=zeros(Nmax,1);
        Ccnt(C_cnt(:,1))=C_cnt(:,k+1);
        %assign to new nodes
        %colonized
        Ccnt(newnodes)=C0;
        %susceptible
        Scnt(newnodes)=1-Ccnt(newnodes);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Scnt,Ccnt,newwardcontamcnt]=agentbasedprob(nl,part,Scnt,Ccnt,para,nodes,appearnodes,alpha,pk,wardcnt,wardcap,wardcontamcnt);
        wardcontam(:,k)=newwardcontamcnt;%update environmental contamination rate
        tempS_cnt(:,k+1)=Scnt(appearnodes);
        tempC_cnt(:,k+1)=Ccnt(appearnodes);
    end
    S_cnt=tempS_cnt;
    C_cnt=tempC_cnt;
end