addpath(genpath(pwd))
format long
a=16384;
nj=a;
ms=a;
xj=sin(-pi*[1:nj]/nj)/2*nj;
cj=exp(1i*pi*[1:nj]/nj)';
k=0:(ms-1);
k=k';

k1=48*cos(pi*[1:ms]/ms);
tol=1e-12;
[fhat1,ffun,U,V]=DeCom_NUFFT1D_III(cj,xj/nj,k1,tol);
fhat2=nufft1d3(nj,xj/nj*2*pi,cj,-1,tol,nk,k1);
error=norm(fhat1-fhat2)/norm(fhat1)
fftconst = -1*1i/ms*2*pi;
fun = @(k,x)exp(fftconst*(k*x'-round(k)*round(x)'));
%K=13;
nufftfun = nufftIII(k1,xj,-1,ms,8*13,tol);
fhat3=nufftfun(cj);
%[U,V] = lowrank(k1,xj,fun,tol,8*13,8*13);
xsub = mod(round(xj),ms)+1;
ksub = mod(round(k1),ms)+1;
r=size(V,2)
Id = sparse(xsub,1:ms,ones(1,ms),ms,ms);
MMM=V.*repmat(cj,[1,r]);
MM=Id*(MMM);
M=fft(MM, [], 1);
fhat =  sum(U.*M(ksub,:),2);
%R=conj(V);
%r=size(V,2);
error1=norm(fhat-fhat2)/norm(fhat)

fid=fopen('./nufftQY/Ur3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',U.');
fclose(fid);
fid=fopen('./nufftQY/Vr3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',V.');
fclose(fid);
fid=fopen('./nufftQY/Ui3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*U.');
fclose(fid);
fid=fopen('./nufftQY/Vi3.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*V.');
fclose(fid);
