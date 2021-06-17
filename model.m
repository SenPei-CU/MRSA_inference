function [incidence,colonized,intro,trans,envir,comm,introc,transc,envirc,rec_inf,rec_col,rec_introc,rec_envirc,appear,newinfection,wardcontam]=model(x,rec_inf,rec_col,rec_introc,rec_envirc,appear,t,tmstep,pararange,para_ens,wardcontam,tstart)
%Run agent-based models for one week
%x:beta - transmission rate; I0 - new introductory infection; C0=gamma - new
%introductory colonization; theta - spillover rate; D - environmental clearance period
%%%%%%input
%(note infected nodes are observed)
%rec_inf - infected nodes at t-1; 
%rec_col - colonized nodes at t-1
%rec_introc - introduced colonization before t
%appaer - nodes ever appeared at t-1
%t - start week
%tmstep - simulation time (week)
%pararange - parameter range
%para_ens - disease parameters for each individual drawn from distributions
%wardcontam - the environmental force of infection
%tstart - start time of simulation
%%%%%%output
%incidence - weekly incidence;(intro+comm+trans)
%colonized - weekly colonization;(introc+transc)
%intro - weekly introduced infection;
%trans - weekly infection of people colonzied in hospitals; 
%comm - weekly infection of people colonized outside hospitals;
%introc - weekly introduced colonization; 
%transc - weekly transmitted colonization
%rec_inf - infected nodes at t; 
%rec_col - colonized nodes at t
%rec_introc - introduced colonization before t+1
%appear - nodes ever appeared at t
%newinfection - newly infected nodes
%wardcontam - updated environmental force of infection

num_ens=size(x,2);
%transform D to lambda
x(end,:)=ones(1,num_ens)./x(end,:);
Nmax=743599;
num_times=1;
incidence=zeros(num_times,num_ens);
colonized=zeros(num_times,num_ens);
intro=zeros(num_times,num_ens);
trans=zeros(num_times,num_ens);
envir=zeros(num_times,num_ens);
comm=zeros(num_times,num_ens);
introc=zeros(num_times,num_ens);
transc=zeros(num_times,num_ens);
envirc=zeros(num_times,num_ens);

%range for C
Cl=pararange(4,1);Cu=pararange(4,2);
load inward
load wardcap
for t1=t:t+tmstep-1
    %load network structure
    filePattern = fullfile(pwd,'net/net');
    load([filePattern,num2str(t1),'.mat']);
    if t1==tstart
        appear(unique(nl(nl(:,2)>0,1)))=1;
    end
    nodes=find(deg>0);
    wardcnt=inward{t1};
    oldnodes=find(appear>0);
    appear(deg>0)=1;
    appearnodes=find(appear>0);
    num_existingnodes=length(appearnodes);
    newnodes=setdiff(appearnodes,oldnodes);
    num_newnodes=length(newnodes);
    for k=1:num_ens
        I0=x(2,k);%importation rate of infected patients, set as 0
        C0=x(3,k);%importation rate of colonized patients
        para=x(1:5,k);%beta0,I0,C0,theta,lambda
        wardcontamcnt=wardcontam(:,k);
        alpha=para_ens(:,1);%decolonization rate
        pk=para_ens(:,2);%observation rate
        mu=para_ens(:,3);%recovery rate after observed, set as alpha
        statecnt=zeros(Nmax,1);
        col_intro=zeros(Nmax,1);
        col_envir=zeros(Nmax,1);
        %recover lastest state
        statecnt(rec_inf{1,k})=2;%assign infection
        statecnt(rec_col{1,k})=1;%assign colonization
        col_intro(rec_introc{1,k})=1;%assign colonization outside hospitals
        col_envir(rec_envirc{1,k})=1;%assign colonization due to environmenal contamination
        
        if t1==tstart%the first step, assign colonized people
            C=Cl+rand(num_existingnodes,1)*(Cu-Cl);
            v=rand(num_existingnodes,1);
            statecnt(appearnodes(v<C))=1;%assign colonized people
        end
        %assign to new nodes
        %colonized
        v=rand(num_newnodes,1);
        statecnt(newnodes(v<C0))=1;
        col_intro(newnodes(v<C0))=1;%record colonization outside hospitals
        colonized(1,k)=colonized(1,k)+sum(v<C0);%colonized
        introc(1,k)=introc(1,k)+sum(v<C0);%introduced colonization
        %infected
        v=rand(num_newnodes,1);
        statecnt(newnodes(v<I0))=2;
        %new infections due to introduction
        incidence(1,k)=incidence(1,k)+sum(v<I0);
        intro(1,k)=intro(1,k)+sum(v<I0);%introducted infection
        newinfintro=newnodes(v<I0);
        [newstate,stat,newinf,newwardcontamcnt,col_envir]=agentbasedsimulation(nl,part,statecnt,para,nodes,appearnodes,col_intro,col_envir,alpha,pk,mu,wardcnt,wardcap,wardcontamcnt);
        wardcontam(:,k)=newwardcontamcnt;%update environmental contamination rate
        
        %stat: 1,trans colonized; 2,envir colonized; 3,infection from intro col;
        %4,infection from trans col; 5,infection from envir col
        
        %find new infections due to transmission and colonization outside
        %hospitals
        trans(1,k)=trans(1,k)+stat(4);%transmitted infection
        comm(1,k)=comm(1,k)+stat(3);
        envir(1,k)=envir(1,k)+stat(5);
        incidence(1,k)=incidence(1,k)+stat(3)+stat(4)+stat(5);
        %find new colonization due to transmission
        transc(1,k)=transc(1,k)+stat(1);%transmitted colonization
        colonized(1,k)=colonized(1,k)+stat(1);
        %find new colonization due to environmental contamination
        envirc(1,k)=envirc(1,k)+stat(2);
        colonized(1,k)=colonized(1,k)+stat(2);
        %update col_introc
        col_intro(newstate~=1)=0;
        %record new state in next time step
        rec_inf{1,k}=find(newstate==2);
        rec_col{1,k}=find(newstate==1);
        rec_introc{1,k}=find(col_intro==1);
        rec_envirc{1,k}=find(col_envir==1);
        %new infection
        newinfection=find(newinf==1);
        newinfection=union(newinfection,newinfintro);
    end
end
