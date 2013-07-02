function [m_pi,m_rho,g_pi,f_piA,m_pcac,tgOm,omega,ZAovZV,ZV,ZA,f_piVps,f_piV] ...
    = Qfit(fullcorr,L,tmin,fitonlast,ntc,mu,sm,sameend,factor)
eval(['addpath /Users/luigi/Work/lib/Matlab/']);
% fullcorr: row (or matrix) of correlators. Indices are expected to have te following order (from slow to fast):
% time(#=T), corr_type(#=ntc), smearing_type (#=sm), conf_index (#=nconf)              [gwc style]
% tmin: minimal time to start cosh/sinh/tanh fits
% fitonlast: range to fit ratios.
% mu: twisted mass.
% L=[Lx, Ly, Lz, T]. 
% sameend: ism index of the correlator which are both smeared or 
% both non-smeared (ll if available, otherwise ss). Usually sameend=1;
% factor (typically the correlator need to be rescaled by a factor = L^3, 2*kappa^2, L^3*2*kappa^2, whatever)

global m_guess xi2_guess rap_guess count1 count2

% usual space-time stuff
D=length(L);
V=prod(L);
sV=prod(L(1:D-1));
T=L(D);
NS=2^(D/2);

count1=0;count2=0;
% doc:do coshfit. dot: do tansh-fit. dor: do ratios. doplot: do plot. (1=yes, 0=no) 
doc=1;dot=1;dor=1;doplot=0;check_therm=0;check_eff=0;

% Type of correlators present inside fullcorr (the order has to match extract_corr)
% iX is the index sX is the signature (cosh/sinh) of the corr.
iPP=1; iAP=2; iPA=3; iAA0=4; iVV0=5; iPV=6; iVP=7; iAV0=8; iVA0=9; iSS=10;
sPP=0; sAP=1; sPA=1; sAA0=0; sVV0=0; sPV=1; sVP=1; sAV0=1; sVA0=1; sSS=0;
i44a=11;iVVs=12;iAAs=13;i4Va=14;iV4a=15;i4Aa=16;iA4a=17;iVAs=18;iAVs=19;iBB=20;
s44a=0; sVVs=0; sAAs=0; s4Va=1; sV4a=1; s4Aa=1; sA4a=1; sVAs=1; sAVs=1; sBB=0;

% derived correlators
iPPmSS=21;iAmVAmV=22;idAP=23;idPA=24;
sPPmSS=0;sAmVAmV=0;sdAP=0;sdPA=0;
delta_ntc=4;

% tell what is missing
iVpsP=[]; iVpsA=[]; iVpsV0=[];

% sign of cosh/sinh
fullS = [sPP sAP sPA sAA0 sVV0 sPV sVP sAV0 sVA0 sSS   s44a sVVs sAAs s4Va sV4a s4Aa sA4a sVAs sAVs sBB ...
	sPPmSS sAmVAmV sdAP sdPA];

% correlator to use:
ctocf_pi = [iPP,iAP,iPA,iAA0];
ctocf_rho= [iVVs]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTINGS ONLY ABOVE
ctocf= [ctocf_pi,ctocf_rho];
lcp=length(ctocf_pi);
lcr=length(ctocf_rho);
nconf=prod(size(fullcorr))/(T*ntc*sm);
%  order needed here:               [nconf T sm ntc]
fulldata = permute(reshape(fullcorr,T,ntc,sm,nconf),[4 1 3 2]); % gwc style
% fulldata = permute(reshape(fullcorr,T,nconf,sm,ntc),[2 1 3 4]); % nf3qcd style

%  rescaling? (typically it was: factor = L^3 ,  2*kappa^2, L^3*2*kappa^2...)
fulldata = fulldata*factor;

% add the derived correlators:
fulldata(:,:,:,iPPmSS)=fulldata(:,:,:,iPP)-fulldata(:,:,:,iSS);
fulldata(:,:,:,iAmVAmV)=fulldata(:,:,:,iAA0)+fulldata(:,:,:,iVV0)-fulldata(:,:,:,iAV0)-fulldata(:,:,:,iVA0);
fulldata(:,1:T-1,:,idAP)=diff(fulldata(:,:,:,iAP),1,2); % 1st derivative along the 2nd index (=T) 
fulldata(:,T,:,idAP)=fulldata(:,1,:,iAP)-fulldata(:,T,:,iAP);
fulldata(:,1:T-1,:,idPA)=diff(fulldata(:,:,:,iPA),1,2); % 1st derivative along the 2nd index (=T) 
fulldata(:,T,:,idPA)=fulldata(:,1,:,iPA)-fulldata(:,T,:,iPA);
ntc=ntc+delta_ntc;

% folding
for itc=1:ntc
  folded(:,[1:T/2-1],:,itc)= ...
      0.5*(fulldata(:,[1:T/2-1]+1,:,itc)+ ((-1)^fullS(itc))*fulldata(:,[T-1:-1:T/2+1]+1,:,itc));
  folded2(:,[1:T-1]+1,:,itc)= ...
      0.5*(fulldata(:,[1:T-1]+1,:,itc)+ ((-1)^fullS(itc))*fulldata(:,[T-1:-1:1]+1,:,itc));
end

t=[tmin:T-tmin+1];
%t=[tmin:T/2];
tf=[tmin-1:T/2-1];

if (check_therm==1)
  figure(10)
  hold on
  jj=0;
  for jt=[T/2-5:T/2-3]
    for jc=[1:3]
      jj=jj+1;
      subplot(length([T/2-5:T/2-3]),length([1:3]),jj)
      plot(fulldata(:,jt,ls,jc))
    end
  end
end
%%%%% Compute Coshfit 
if(doc==1)
  for ism= 1:sm;
    itcj=0;
    for itc= ctocf_pi;
      itcj=itcj+1;
      m_guess=5e-1;
      xi2_guess=2*m_guess*exp(m_guess*(tmin-1))*mean(folded2(:,tmin,ism,itc),1);
      data = squeeze(folded2(:,t,ism,itc));
      for v=t;
	[value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,v-tmin+1);
	co(v,ism,itc)=value;
	dco(v,ism,itc)=dvalue;
      end 
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@coshfit,dco(t,ism,itc)',1,fullS(itc),t,T);
      ampl(ism,itc)=value;
      dampl(ism,itc)=dvalue;
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@coshfit,dco(t,ism,itc)',2,fullS(itc),t,T);
      mass(ism,itc)=value;
      dmass(ism,itc)=dvalue;
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@coshfit,dco(t,ism,itc)',3,fullS(itc),t,T);
      sqrtampl(ism,itc)=value;
      dsqrtampl(ism,itc)=dvalue;
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@coshfit,dco(t,ism,itc)',4,fullS(itc),t,T);
      sqrtamplovm(ism,itc)=value;
      dsqrtamplovm(ism,itc)=dvalue;
      %  tr=[tmin:T-tmin-1];
      tr=[tmin:T/2-1];
      meff(tr)=log(co(tr,ism,itc)./co(tr+1,ism,itc));
      dmeff(tr)=sqrt( (dco(tr,ism,itc)./co(tr,ism,itc)).^2 + (dco(tr+1,ism,itc)./co(tr+1,ism,itc)).^2   );
      %%%%%%%%%% Plot Coshfit
      if(doplot==1)
	figure(1)
	subplot(sm,lcp,(ism-1)*lcp+itcj)
	hold on
	plot(t,coshfun(t,T,ampl(ism,itc),mass(ism,itc),fullS(itc)),'b-')
	errorbar(t,co(t,ism,itc),dco(t,ism,itc),'rx')
	hold off
      end
      if(check_eff==1)
	figure(100)
	subplot(sm,lcp,(ism-1)*lcp+itcj)
	hold on
	errorbar(tr,meff(tr),dmeff(tr))
      end
    end
  end
  for ism= 1:sm;
    itcj=0;
    for itc= ctocf_rho;
      itcj=itcj+1;
      m_guess=7e-1;
      xi2_guess=2*m_guess*exp(m_guess*(tmin-1))*mean(folded2(:,tmin,ism,itc),1);
      data = squeeze(folded2(:,t,ism,itc));
      for v=t;
	[value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,v-tmin+1);
	co(v,ism,itc)=value;
	dco(v,ism,itc)=dvalue;
      end 
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@coshfit,dco(t,ism,itc)',1,fullS(itc),t,T);
      ampl(ism,itc)=value;
      dampl(ism,itc)=dvalue;
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@coshfit,dco(t,ism,itc)',2,fullS(itc),t,T);
      mass(ism,itc)=value;
      dmass(ism,itc)=dvalue;
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@coshfit,dco(t,ism,itc)',3,fullS(itc),t,T);
      sqrtampl(ism,itc)=value;
      dsqrtampl(ism,itc)=dvalue;
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@coshfit,dco(t,ism,itc)',4,fullS(itc),t,T);
      sqrtamplovm(ism,itc)=value;
      dsqrtamplovm(ism,itc)=dvalue;
      %  tr=[tmin:T-tmin-1];
      %tr=[tmin:T/2-1];
      %meff(tr)=log(co(tr,ism,itc)./co(tr+1,ism,itc));
      
      %%%%%%%%%% Plot Coshfit
      if(doplot==1)
	figure(2)
	subplot(sm,lcr,(ism-1)*lcr+itcj)
	hold on
	plot(t,coshfun(t,T,ampl(ism,itc),mass(ism,itc),fullS(itc)),'b-')
	errorbar(t,co(t,ism,itc),dco(t,ism,itc),'rx')
	hold off
	%figure(2)
	%plot(tr,meff(tr))
      end
    end
  end
end
%%%%% Compute Tanhfit 
if(dot==1)
  for ism= 1:sm;
    [ism iAP iPP];
    for v=t;
      data = squeeze(folded2(:,v,ism,[iAP,iPP]));
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@ratio);
      rco(v,ism)=value;
      drco(v,ism)=dvalue;
    end 
    m_guess=5e-1;
    rap_guess=(mean(folded2(:,tmin,ism,iAP),1)/mean(folded2(:,tmin,ism,iPP),1))/...
	tanh(m_guess*(T/2+1-tmin ));
    data = squeeze(folded2(:,t,ism,[iAP,iPP]));
    data=reshape(data,nconf,length(t)*2);
    [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drco(t,ism)',1,t,T);
    rampl(ism)=value;
    drampl(ism)=dvalue;
    [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drco(t,ism)',2,t,T);
    rmass(ism)=value;
    drmass(ism)=dvalue;
    [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drco(t,ism)',3,t,T);
    ramplovmass(ism)=value;
    dramplovmass(ism)=dvalue;
    [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drco(t,ism)',4,t,T);
    mpcac(ism)=value;
    dmpcac(ism)=dvalue;

    [ism iVP iPP];
    for v=t;
      data = squeeze(folded2(:,v,ism,[iVP,iPP]));
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@ratio);
      rcoV(v,ism)=value;
      drcoV(v,ism)=dvalue;
    end 
    m_guess=5e-1;
    rap_guess=(mean(folded2(:,tmin,ism,iVP),1)/mean(folded2(:,tmin,ism,iPP),1))/...
	      tanh(m_guess*(T/2+1-tmin ));
    data = squeeze(folded2(:,t,ism,[iVP,iPP]));
    data=reshape(data,nconf,length(t)*2);
    [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drcoV(t,ism)',1,t,T);
    ramplV(ism)=value;
    dramplV(ism)=dvalue;
    [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drcoV(t,ism)',2,t,T);
    rmassV(ism)=value;
    drmassV(ism)=dvalue;
    [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drcoV(t,ism)',3,t,T);
    ramplovmassV(ism)=value;
    dramplovmassV(ism)=dvalue;
    [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drcoV(t,ism)',4,t,T);
    invZV_indovmu(ism)=value;
    dinvZV_indovmu(ism)=dvalue;
    if(~isempty(iVpsP))
      [ism iVpsP iPP];
      for v=t;
	data = squeeze(folded2(:,v,ism,[iVpsP,iPP]));
	[value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@ratio);
	rcoVps(v,ism)=value;
	drcoVps(v,ism)=dvalue;
      end 
      m_guess=5e-1;
      rap_guess=(mean(folded2(:,tmin,ism,iVpsP),1)/mean(folded2(:,tmin,ism,iPP),1))/...
		tanh(m_guess*(T/2+1-tmin ));
      data = squeeze(folded2(:,t,ism,[iVpsP,iPP]));
      data=reshape(data,nconf,length(t)*2);
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drcoVps(t,ism)',1,t,T);
      ramplVps(ism)=value;
      dramplVps(ism)=dvalue;
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drcoVps(t,ism)',2,t,T);
      rmassVps(ism)=value;
      drmassVps(ism)=dvalue;
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drcoVps(t,ism)',3,t,T);
      ramplovmassVps(ism)=value;
      dramplovmassVps(ism)=dvalue;
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tanhfit,drcoVps(t,ism)',4,t,T);
      invZVps_indovmu(ism)=value;
      dinvZVps_indovmu(ism)=dvalue;
    end
    %%%%%%%%%% Plot Tanhfit
    if(doplot==1)
      figure(3)
      subplot(sm,1,ism)
      hold on
      plot(t,tanhfun(t,T,rampl(ism),rmass(ism)),'b-')
      errorbar(t,rco(t,ism),drco(t,ism),'rx')
      hold off
      figure(4)
      subplot(sm,1,ism)
      hold on
      plot(t,tanhfun(t,T,ramplV(ism),rmassV(ism)),'b-')
      errorbar(t,rcoV(t,ism),drcoV(t,ism),'rx')
      hold off
      if(~isempty(iVpsP))
	figure(5)
	subplot(sm,1,ism)
	hold on
	plot(t,tanhfun(t,T,ramplVps(ism),rmassVps(ism)),'b-')
	errorbar(t,rcoVps(t,ism),drcoVps(t,ism),'rx')
	hold off
      end
    end
  end
end
if(dor==1)
  %%%%%%%%%% Compute Ratios, Omegas and Z's
  for ism =1:sm
    % 1. tan(omega_V) = VP0/AP0
    ob=1;
    for v=tf;
      data = squeeze(folded(:,v,ism,[iVP,iAP]));
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@ratio);
      rat(v,ism,ob)=value; %was: v+tmin-1
      drat(v,ism,ob)=dvalue;
    end 
    fit(ism,ob)=mean(rat(fitonlast,ism,ob),1);
    dfit(ism,ob)=mean(drat(fitonlast,ism,ob),1);
    
    % 2. tan(omega_A) = (AVs - tgOmV_ss AAs)/(VVs + tgOmV_ss AVs)
    ob=2;
    for v=tf;
      data = squeeze(folded(:,v,ism,[iVVs,iAVs,iVAs,iAAs]));
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@tgOmA,rat(v,1,1));
      rat(v,ism,ob)=value;
      drat(v,ism,ob)=dvalue;
    end 
    fit(ism,ob)=mean(rat(fitonlast,ism,ob),1);
    dfit(ism,ob)=mean(drat(fitonlast,ism,ob),1);
    
    % 3. tan(omega) 
    ob=3;
    for v=tf;
      data1 = squeeze(folded(:,v,ism,[iVVs,iAVs,iVAs,iAAs]));
      data2 = squeeze(folded(:,v,sameend,[iVP,iAP]));
      data = [data1,data2];
      docount =sum(fitonlast==v);
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@Gettgomega,docount);
      rat(v,ism,ob)=value;
      drat(v,ism,ob)=dvalue;
    end 
    fit(ism,ob)=mean(rat(fitonlast,ism,ob),1);
    dfit(ism,ob)=mean(drat(fitonlast,ism,ob),1);
    
    % 4. omega 
    ob=4;
    for v=tf;
      data1 = squeeze(folded(:,v,ism,[iVVs,iAVs,iVAs,iAAs]));
      data2 = squeeze(folded(:,v,sameend,[iVP,iAP]));
      data = [data1,data2];
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@Getomega);
      rat(v,ism,ob)=value;
      drat(v,ism,ob)=dvalue;
    end 
    fit(ism,ob)=mean(rat(fitonlast,ism,ob),1);
    dfit(ism,ob)=mean(drat(fitonlast,ism,ob),1);
    
    if(~isempty(iVpsP))
      % 6. Z_V 
      ob=6;
      for v=tf;
	data = squeeze(folded(:,v,ism,[iVpsP,iVP]));
	[value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@ratio);
	rat(v,ism,ob)=value;
	drat(v,ism,ob)=dvalue;
      end 
      fit(ism,ob)=mean(rat(fitonlast,ism,ob),1);
      dfit(ism,ob)=mean(drat(fitonlast,ism,ob),1);
    end
    
    % 7. ZAovZV
    ob=7;
    for v=tf;
      data1 = squeeze(folded(:,v,ism,[iVVs,iAVs,iVAs,iAAs]));
      data2 = squeeze(folded(:,v,sameend,[iVP,iAP]));
      data = [data1,data2];
      [value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@GetZAovZV);
      rat(v,ism,ob)=value;
      drat(v,ism,ob)=dvalue;
    end 
    fit(ism,ob)=mean(rat(fitonlast,ism,ob),1);
    dfit(ism,ob)=mean(drat(fitonlast,ism,ob),1);
    
    if(~isempty(idAP))
      % 5. m_PCAC \chi
      ob=5;
      for v=tf;
	data = squeeze(folded(:,v,ism,[idAP,iPP]));
	[value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@ratio);
	rat(v,ism,ob)=value;
	drat(v,ism,ob)=dvalue;
      end 
      fit(ism,ob)=mean(rat(fitonlast,ism,ob),1);
      dfit(ism,ob)=mean(drat(fitonlast,ism,ob),1);
    end
    if(~isempty(idPA))
      % 5. m_PCAC \chi
      ob=8;
      for v=tf;
	data = squeeze(folded(:,v,ism,[idPA,iPP]));
	[value,dvalue,ddvalue,tauint,dtauint,Qval] = UWerr(data,[],[],0,@ratio);
	rat(v,ism,ob)=value;
	drat(v,ism,ob)=dvalue;
      end 
      fit(ism,ob)=mean(rat(fitonlast,ism,ob),1);
      dfit(ism,ob)=mean(drat(fitonlast,ism,ob),1);
    end
  end
  %%%%%%%%%% Plot Ratios, Omegas and Z's
  if(doplot==1)
    figure(6)
    j=0;
    for ism=1:sm
      for ob=1:4
	j=j+1;
	subplot(sm,4,j)
	hold on
	errorbar(tf,rat(tf,ism,ob),drat(tf,ism,ob),'rx')
	plot(fitonlast,fit(ism,ob)+dfit(ism,ob),'b+')
	plot(fitonlast,fit(ism,ob)-dfit(ism,ob),'b+')
	hold off
	end
    end
    figure(8)
    j=0;
    for ism=1:sm
      for ob=7
	j=j+1;
	subplot(sm,2,j)
	hold on
	errorbar(tf,rat(tf,ism,ob),drat(tf,ism,ob),'rx')
	plot(fitonlast,fit(ism,ob)+dfit(ism,ob),'b+')
	plot(fitonlast,fit(ism,ob)-dfit(ism,ob),'b+')
	hold off
	end
    end
    if(~isempty(idAP))
      figure(7)
      j=0;
      for ism=1:sm
	for ob=5
	  j=j+1;
	  subplot(sm,1,j)
	  hold on
	  errorbar(tf,rat(tf,ism,ob),drat(tf,ism,ob),'rx')
	  plot(fitonlast,fit(ism,ob)+dfit(ism,ob),'b+')
	  plot(fitonlast,fit(ism,ob)-dfit(ism,ob),'b+')
	  hold off
	end
      end
    end
    if(~isempty(idPA))
      figure(17)
      j=0;
      for ism=1:sm
	for ob=8
	  j=j+1;
	  subplot(sm,1,j)
	  hold on
	  errorbar(tf,rat(tf,ism,ob),drat(tf,ism,ob),'rx')
	  plot(fitonlast,fit(ism,ob)+dfit(ism,ob),'b+')
	  plot(fitonlast,fit(ism,ob)-dfit(ism,ob),'b+')
	  hold off
	end
      end
    end
    if(~isempty(iVpsP))
      figure(9)
      j=0;
      for ism=1:sm
	for ob=6
	  j=j+1;
	  subplot(sm,1,j)
	  hold on
	  errorbar(tf,rat(tf,ism,ob),drat(tf,ism,ob),'rx')
	  plot(fitonlast,fit(ism,ob)+dfit(ism,ob),'b+')
	  plot(fitonlast,fit(ism,ob)-dfit(ism,ob),'b+')
	  hold off
	end
      end
    end
  end
end
%%%%%%%%%%%%%%   Output Results
count1;
count2;
m_pi=[];m_rho=[];m_pcac=[];g_pi=[];f_piA=[];f_piVps=[];f_piV=[];tgOm=[];omega=[];ZAovZV=[];ZV=[];ZA=[];
for ism=1:sm
  % m_pi: PP    AP   PA   AA   AP/PP
  m_pi(:,:,ism)=[ mass(ism,iPP), mass(ism,iAP), mass(ism,iPA), mass(ism,iAA0), rmass(ism);...
        	     dmass(ism,iPP),dmass(ism,iAP),dmass(ism,iPA), dmass(ism,iAA0),drmass(ism)];
  % m_rho: VVs
  m_rho(:,:,ism)=[mass(ism,iVVs); ...
              dmass(ism,iVVs)];
end
% choose m_pi 
I = [1 2];
tmpi=m_pi(:,I,:);
tmpi=reshape(tmpi,2,size(tmpi,2)*sm);
tmrho=reshape(m_rho,2,size(m_rho,2)*sm);
avm_pi(1,1)=sum(tmpi(1,:)./(tmpi(2,:)).^2)/sum(1./(tmpi(2,:)).^2);
avm_pi(2,1)=1/sqrt(sum(1./(tmpi(2,:)).^2));
avm_rho(1,1)=sum(tmrho(1,:)./(tmrho(2,:)).^2)/sum(1./(tmrho(2,:)).^2);
avtmrho(2,1)=1/sqrt(sum(1./(tmrho(2,:)).^2));
clear tmpi tmrho
for ism=1:sm
  % m_pcac: rAP dAP/PP dPA/PP
  m_pcac(:,:,ism)=[mpcac(ism),0.5*fit(ism,5),0.5*fit(ism,8);...
		   dmpcac(ism),0.5*dfit(ism,5),0.5*dfit(ism,8)];
end
for ism=1:sm
  % g_pi: PP
  if (ism==sameend)
    g_pi(:,1,ism) = [sqrtampl(ism,iPP);...
		     dsqrtampl(ism,iPP)];
  else
    g_pi(1,1,ism) = ampl(ism,iPP)/sqrtampl(sameend,iPP);
    g_pi(2,1,ism) = abs(temp(1,1))*sqrt( (dampl(ism,iPP)/ampl(ism,iPP))^2 + (dsqrtampl(sameend,iPP)/sqrtampl(sameend,iPP))^2 );
  end
end
for ism=1:sm
  % f_pi_Chi_A: AA rAP
  if (ism==sameend)
    temp=[sqrtamplovm(ism,iAA0);...
	  dsqrtamplovm(ism,iAA0)];
  else
    temp(1,1)=ampl(ism,iAA0)/(sqrtampl(sameend,iAA0)*avm_pi(1,1));
    temp(2,1)=abs( temp(1,1) )*...
	      sqrt((dampl(ism,iAA0)/ampl(ism,iAA0))^2+...
		   (dsqrtampl(sameend,iAA0)/sqrtampl(sameend,iAA0))^2+...
		   (avm_pi(2,1)/avm_pi(1,1))^2);
  end
  temp1=[ramplovmass(ism)*g_pi(1,1,ism);...
	 sqrt((dramplovmass(ism)*g_pi(1,1,ism))^2+(ramplovmass(ism)*g_pi(2,1,ism))^2)];
  f_piA(:,:,ism)=[temp,temp1];
  clear temp temp1

  % f_pi_Chi_Vps: (from <|Vcons|pi>/m_pi ) rVpsPss rVpsPls
  if(~isempty(iVpsP))
    f_piVps(:,:,ism)=[ramplovmassVps(ism)*g_pi(1,1,ism);...
		      sqrt((dramplovmassVps(ism)*g_pi(1,1,ism))^2+(ramplovmassVps(ism)*g_pi(2,1,ism))^2)];
  else
    f_piVps(:,:,ism)=nan(2,1);
  end
  
  % f_pi_Chi_V: (from <|V|pi>/m_pi )   rVPss rVPls
  temp=[ramplovmassV(ism)*g_pi(1,1,ism);...
	 sqrt((dramplovmassV(ism)*g_pi(1,1,ism))^2+(ramplovmassV(ism)*g_pi(2,1,ism))^2)];
  temp1(1,1)=mu*g_pi(1,1,ism)/avm_pi(1,1)^2;
  temp1(2,1)=abs(temp1(1,1)) * sqrt((g_pi(2,1,ism)/g_pi(1,1,ism))^2+(2*avm_pi(2,1)/avm_pi(1,1))^2);

  f_piV=[temp,temp1];
  clear temp temp1
end
for ism=1:sm
  %tgO:  tgomega_ss tgomega_ls tgOVss tgOVls tgOAss tgOAls
  tgOm(:,:,ism) = [  fit(ism,3), fit(ism,1), fit(ism,2);...
		     dfit(ism,3),dfit(ism,1),dfit(ism,2)];
  % omega_ss/pi omega_ls/pi
  omega(:,:,ism) = [fit(ism,4);dfit(ism,4)]./pi;
  %ZAovZV_ss ZAovZV_ls 
  ZAovZV(:,:,ism)= [fit(ism,7);dfit(ism,7)];
  if(~isempty(iVpsP))
    %ZV     ZV_c_ss ZV_c_ls  ZV_ind_ss ZV_ind_ls
    ZV(:,:,ism) = [fit(ism,6),...
		   mu/invZV_indovmu(ism);...
		   dfit(ism,6),... 
		   mu*dinvZV_indovmu(ism)/(invZV_indovmu(ism))^2];
  else
    %ZV  ZV_ind_ss ZV_ind_ls
    ZV(:,:,ism) = [mu/invZV_indovmu(ism);...
		   mu*dinvZV_indovmu(ism)/(invZV_indovmu(ism))^2];
  end
end
for ism=1:sm
  if(~isempty(iVpsP))
    %ZA  ZA_mpcac_ss ZA_mpcac_ls ZA_om_c_ss ZA_om_c_ls ZA_om_ind_ss ZA_om_ind_ls    
    %  ... later also ZA_fpi]
    ZA(:,:,ism) = [mu/(tgOm(1,1,ism)*m_pcac(1,1,ism)), ...
		 ZAovZV(1,1,ism)*ZV(1,1,ism), ...
		 ZAovZV(1,1,ism)*ZV(1,2,ism); ...
		 abs(mu/(tgOm(1,1,ism)*m_pcac(1,1,ism)))*sqrt((tgOm(2,1,ism)/tgOm(1,1,ism))^2+(m_pcac(2,1,ism)/m_pcac(1,1,ism))^2),...
		 sqrt( (ZAovZV(2,1,ism)*ZV(1,1,ism))^2 +(ZAovZV(1,1,ism)*ZV(2,1,ism))^2),...
		 sqrt( (ZAovZV(2,1,ism)*ZV(1,2,ism))^2 +(ZAovZV(1,1,ism)*ZV(2,1,ism))^2)];
    else
    ZA(:,:,ism) = [mu/(tgOm(1,1,ism)*m_pcac(1,1,ism)), ...
		 ZAovZV(1,1,ism)*ZV(1,1,ism); ...
		 abs(mu/(tgOm(1,1,ism)*m_pcac(1,1,ism)))*sqrt((tgOm(2,1,ism)/tgOm(1,1,ism))^2+(m_pcac(2,1,ism)/m_pcac(1,1,ism))^2),...
		 sqrt( (ZAovZV(2,1,ism)*ZV(1,1,ism))^2 +(ZAovZV(1,1,ism)*ZV(2,1,ism))^2)];
  end      
end % sm 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS
%%%%%%%%% cosh
function [res f] = coshfit(corr,dcorr,P,s,t,T)
global m_guess xi2_guess rap_guess
start_guess = [xi2_guess,m_guess];
model=@(par)coshdev(par,corr,dcorr,s,t,T);
opt=optimset('MaxFunEvals',1000000,'MaxIter',1000000); % default 1000
[resi,f,ef]=fminsearch(model,start_guess); % opt
if((P==1))
  res=resi(1);
elseif(P==2)
  res=abs(resi(2));
elseif(P==3)
  res=sqrt(abs(resi(1)));
elseif(P==4)
  res=sqrt(abs(resi(1)))/abs(resi(2));
end
if (ef<1)
  [ef,f,1]
  res=nan;
  f=nan;
end
function sse = coshdev(par,corr,dcorr,s,t,T)
xi=par(1);
m=par(2);
thcorr =  coshfun(t,T,xi,m,s); 
sse= sum(((thcorr - corr)./(dcorr+1e-10)).^2);
function corr = coshfun(t,T,xiq,m,s)
m=abs(m);
corr =  xiq/(2*m) * (exp(-m*(t-1)) + (-1)^s * exp(-m*(T-(t-1))));
%%%%%%%% tanh
function [res f] = tanhfit(corr,dcorr,P,t,T)
global m_guess xi2_guess rap_guess
start_guess = [rap_guess,m_guess];
model=@(par)tanhdev(par,corr,dcorr,t,T);
opt=optimset('MaxFunEvals',1000000,'MaxIter',1000000); % default 1000
[resi,f,ef]=fminsearch(model,start_guess);% opt
if((P==1))
  res=resi(1);
elseif(P==2)
  res=abs(resi(2));
elseif(P==3)
  res=resi(1)/abs(resi(2));
elseif(P==4)
  res=0.5*resi(1)*abs(resi(2));
end
if (ef<1)
  [ef,f,2]
  res=nan;
  f=nan;
end
function sse = tanhdev(par,corr,dcorr,t,T)
rap=par(1);
m=par(2);
thcorr =  tanhfun(t,T,rap,m);
lt=length(t); 
sse= sum(((thcorr - (corr(1:lt)./corr(lt+1:2*lt) ) )./(dcorr(1:lt)+1e-10)).^2);
function ratcorr = tanhfun(t,T,rap,m)
m=abs(m);
ratcorr = rap * tanh(m*(T/2-(t-1)));
%%%%%%%% ratios 
function rr = ratio(xx)
rr = xx(1)./xx(2);
%%%%%%%% Omega's
function tg= tgOmA(xx,tgV);
tg = (xx(2) - tgV .* xx(4))./(xx(1) + tgV .* xx(3));
function tg= tgOmSW(xx,tgV);
tg = (xx(3) -  xx(2))./(xx(1) - xx(4));
function omega= Getomega(xx);
tgOV = xx(5)./xx(6);
tgOA = (xx(2) - tgOV .* xx(4))./(xx(1) + tgOV .*xx(3));
sv=sign(tgOV);
sa=sign(tgOA);
omega = atan(sv*sqrt(abs(tgOV*tgOA)));
if(omega<0)
  omega=omega+pi;
end
function tgomega= Gettgomega(xx,count);
global count1 count2
tgOV = xx(5)./xx(6);
tgOA = (xx(2) - tgOV .* xx(4))./(xx(1) + tgOV .*xx(3));
sv=sign(tgOV);
sa=sign(tgOA);
if(count==1)
  count1=count1+1;
  if(sa~=sv)
    count2=count2+1;
  end
end
tgomega = sv*sqrt(abs(tgOV*tgOA));
function zav= GetZAovZV(xx);
tgOV = xx(5)./xx(6);
tgOA = (xx(2) - tgOV .* xx(4))./(xx(1) + tgOV .*xx(3));
zav = sqrt(abs(tgOV./tgOA));
