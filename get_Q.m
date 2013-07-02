eval(['addpath /Users/luigi/Work/_lib/Matlab/']);
% read correlator files produced by get_corr
% as listed in table_run.dat; 
% Comment out lines in table_run.dat to skip some points. 

eval(['load table_run.dat']);
eval(['load table_Q.mine3']);
table_Q1=table_Q;
jj=0;
for n=1:size(table_run,1);
  n
  be=table_run(n,1);
  mu=table_run(n,2);
  L=table_run(n,3);
  T=table_run(n,4);
  ka=table_run(n,5);
  lh=table_run(n,6);
  nc=table_run(n,7);
  ns=table_run(n,8);
  tmin=ceil(T/5);
  fitonlast=[T/2-5:T/2-3];
  if((be==0.74)&(mu==0.0075)&(ka==0.159)&(lh==1))
    fitonlast=[T/2-7:T/2-6];
  end
  if((be==0.67)&(mu==0.01)&(ka==0.167)&(lh==1))
    fitonlast=[T/2-5:T/2-4];
  end
  if((be==5.3)&(mu==0.008)&(ka==0.16735)&(lh==1))
    fitonlast=[T/2-7:T/2-6];
  end
  if((be==5.1)&(mu==0.013)&(ka==0.1758)&(lh==1))
    fitonlast=[T/2-5:T/2-4];
  end

  if(lh==-1)
    eval(['load ../correlators/corr.m' num2str(mu) '_b' num2str(be) '_L' num2str(L) 'T' num2str(T) '_k' num2str(ka) '_low.dat']);
  elseif (lh==1)
    eval(['load ../correlators/corr.m' num2str(mu) '_b' num2str(be) '_L' num2str(L) 'T' num2str(T) '_k' num2str(ka) '_high.dat']);
  end


[m_pi,m_rho,g_pi,f_pi,m_pcac,tgOm,omega,ZAovZV,ZV,ZA,junk1,f_piVps,junk2,f_piV] = Qfit(corr,L,T,tmin,fitonlast,nc,ns,mu,ka);

res(:,:,n)=[m_pi,m_rho,g_pi, f_pi, m_pcac,tgOm, omega,ZAovZV,ZV,ZA,junk1,f_piVps,junk2,f_piV];

ind=find((table_Q1(:,1)==L)&(table_Q1(:,2)==T)&(table_Q1(:,3)==lh)&(table_Q1(:,4)==be)&...
	 (table_Q1(:,5)==ka)&(table_Q1(:,6)==mu));
if (isempty(ind))
  jj=jj+1;
  ind = size(table_Q1,1)+jj;
end
table_Q1(ind,1:6)=[L, T, lh, be, ka, mu];
%%%%%% table_Q1(ind,7:10)=[r0th, dr0th, r0hyp, dr0hyp];
table_Q1(ind,11:2:10+2*50)=res(1,:,n);
table_Q1(ind,12:2:11+2*50)=res(2,:,n);

%%%% Colums in table_Q, res (Qfit()) have the following meaning:
% table_Q,res
% 1:6           L  T  L/H  beta  kappa  mu
% 7:10         (r0):      r0-thin  r0-hyp
% 11:28,  1: 9 (m_pi) :  av PPss   PPls   APss APls AAss AAls AP/PPss AP/PPls  
% 29:34, 10:12 (m_rho):  av VVi_ss VVi_ls
% 35:38, 13:14 (g_pi):   PPss   PPls
% 39:46, 15:18 (f_pi_Chi):   AAss    AAls   rAPss     rAPls
% 47:54, 19:22 (m_pcac_Chi_1,2): rAP_ss  rAP_ls  dAP/PP_ss dA/PP_ls
% 55:66, 23:28 (tgO):   tgomega_ss tgomega_ls tgOVss tgOVls tgOAss tgOAls
% 67:70, 29:30 (omega): omega_ss/pi omega_ls/pi
% 71:74, 31:32 (ZA/ZV): ZAovZV_ss ZAovZV_ls 
% 75:82, 33:36 (ZV):    ZV_c_ss ZV_c_ls  ZV_ind_ss ZV_ind_ls
% 83:94, 37:42 (ZA):    ZA_mpcac_ss ZA_mpcac_ls ZA_om_c_ss ZA_om_c_ls ZA_om_ind_ss ZA_om_ind_ls    
% 95:102,43:46(f_pi_Chi_Vps): VVss VVls rVPss rVPls
% 103:110,47:50(f_pi_Chi_V):  VVss VVls rVPss rVPls

end

save table_Q.temp table_Q1 -ASCII

