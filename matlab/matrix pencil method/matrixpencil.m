function [freq,amp,varargout]=matrixpencil(x,cut_off)
%MATRIXPENCIL  Matrix Pencil Analysis of Signal Spectrum
% [FREQ,AMP]=MATRIXPENCIL(X,CUT_OFF) Returns frequencies and amplitudes 
% of main spectral components of input signal given in vector X.  Complex 
% angular frequencies are returned in vector FREQ and amplitudies in vector
% AMP. Angular frequencies are given in units of 1/dt where dt is the time 
% or space discretization step of signal. Amplitudes are returned in units 
% of signal a CUT_OFF parameter set the cut off value for singular values 
% calculated in script; raising the CUT_OFF increase number of found
% frequencies with small amplitude but may lead to non-convergence of 
% method or wrong results particulary for components with small amplitude. 
% Deacreasing CUT_OFF value eliminates the small amplitiude results from 
% output. Default value of CUT_OFF is 5.
%
% For input signal X and length N output vectors FREQ and AMP 
%
%       N
% X(t)=sum  AMP(n)*exp(j*FREQ(n)*(t-1)) + NOISE(t)
%      n=1
%
% [FREQ,AMP,FLAG,RELRES]=MATRIXPENCIL(...) returns also convergence flag 
% and relative residual norm from least square solution for amplitudes AMP;
% see LSQR for further details
% FREQ=MATRIXPENCIL(...) returns only the frequencies
%                                              v 3.0     W.M.Saj July 2007
% Literature:
% T.K.Sarkar and O. Pereira
% "Using the Matrix Pencil Method to Estimate the Parameters of a Sum 
% of Complex Exponentials" 
% IEEE Antennas and Propagation Magazine, Vol.27, No 1 Feb 1995 

%% checking output parameters number
if (nargout>4) | (nargout<2)
    error('improper number of output parameters');return;
end

%% setting default value for cut_off
if nargin==1
    cut_off=5;
end

N=length(x);
L=ceil(N/3);

%%Hankel data matrix
Y=x(hankel(1:N-L,N-L:N));

%% SVD
[U,S,V] = svds(Y,N-L);

%% FILTERING
%% getting singular values
singval=diag(S);
maxsingval=singval(1);
for indx=2:L
    M=indx;
    if log10(abs(singval(indx)/maxsingval))<-cut_off  break; end
end

Sp(1:M,1:M)=S(1:M,1:M);clear S;
Vp=V(:,1:M);clear V;
Up=U(:,1:M);clear U;

%% Y1 & Y2
[c,r]=size(Vp);
Y1=Up*Sp*(Vp(1:r-1,:)');
Y2=Up*Sp*(Vp(2:r,:)');
clear Up;clear Sp;clear Vp;

%% spectrum of matrix pencil
z=eig(pinv(Y1)*Y2,'balance');
freq=log(z)/(sqrt(-1));freq=conj(freq);

%% residues
A=ones(N,length(z));
for indx=2:N
    A(indx,:)=A(indx-1,:).*z.';
end
[amp,flag,relres]=lsqr(A,x.',[],[],[],[],pinv(A)*(x.'));

if nargout>2
    varargout(1)={flag};
    if nargout>3
        varargout(2)={relres};
    end
end
