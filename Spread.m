classdef Spread
    % SPREAD is a class that calculate spread option price and Greeks by
    % FFT. There is also plot function within this class.

    % User need to download the complex Gamma function from
    % https://www.mathworks.com/matlabcentral/fileexchange/3572-gamma
    % by Paul Godfrey and rename it "cgamma".
    
    % The function Spread.Price and Spread.Greek calculate price and Greek    
    % Parameters - (S1,S2,K,r,sigma1,sigma2,rho,tau)
    % S1 - asset 1 intital price
    % S1 - asset 2 intital price
    % K - strike Price
    % r - risk free rate
    % sigma1 - volatility of asset 1
    % sigma2 - volatility of asset 2
    % rho - correlation between the assets
    % tau - time to maturity
    
    % The function Spread.plot will create a 3-D surface plot the Greeks
    % against the intital price of the two assets
    % Parameters - (K,r,sigma1,sigma2,rho,tau,N)
    % N - plotting range (1 to N)
    
    properties
    end
    
    methods(Static)    
        function [p1,p2,C,HH1,HH2,HT,HHH11,HHH12,HHH22]=FFTMatrix(S1,S2,K,r,sigma1,sigma2,rho,tau)
            N=128;
            u_bar=40;
            e1=-3;
            e2=1;
            
            X1 = log(S1/K);
            X2 = log(S2/K);
            
            l = 0:N-1;
            
            u_bar1=Spread.utest(u_bar,X1,N);
            u_bar2=Spread.utest(u_bar,X2,N);
            
            eta1 = 2*u_bar1 / N;
            eta2 = 2*u_bar2 / N;
            eta_star1 = pi / u_bar1;
            eta_star2 = pi / u_bar2;
            
            u1 = -u_bar1 + eta1*l;
            u2 = -u_bar2 + eta2*l;
            x1 = -0.5*N*eta_star1 + eta_star1*l;
            x2 = -0.5*N*eta_star2 + eta_star2*l;
            
            H=zeros(N,N);
            HT=zeros(N,N);
            HH1=zeros(N,N);
            HH2=zeros(N,N);
            HHH11=zeros(N,N);
            HHH12=zeros(N,N);
            HHH22=zeros(N,N);
            C=zeros(N,N);
            
            rvec = [r,r]';
            sigma = [sigma1^2, sigma2^2]';
            u_vec = rvec-0.5*sigma;
            Sig = [sigma1^2, rho*sigma1*sigma2; rho*sigma1*sigma2, sigma2^2];
            
            for i=0:N-1
                for j=0:N-1
                    u=[u1(i+1)+1i*e1,u2(j+1)+1i*e2];
                    mu=u*u_vec;
                    Sigma=u*Sig*u';
                    phi = Spread.Phi(mu,Sigma,tau);
                    phat = Spread.P_hat(u);
                    H(i+1, j+1)=(-1)^(i+j) * phi * phat;
                    HT(i+1, j+1)=(1i*mu-0.5*Sigma-r)*H(i+1, j+1);
                    HH1(i+1, j+1)=(1i*u1(i+1)-e1)*H(i+1, j+1);
                    HH2(i+1, j+1)=(1i*u2(j+1)-e2)*H(i+1, j+1);
                    HHH11(i+1, j+1)=(1i*u1(i+1)-e1)^2*H(i+1, j+1);
                    HHH12(i+1, j+1)=(1i*u1(i+1)-e1)*(1i*u2(j+1)-e2)*H(i+1, j+1);
                    HHH22(i+1, j+1)=(1i*u2(j+1)-e2)^2*H(i+1, j+1);
                    
                    C(i+1, j+1)=(-1)^(i+j)*exp(-e1*x1(i+1)-e2*x2(j+1))*eta1*eta2*(N/(2*pi))^2;
                end
            end
            z1 = abs((X1+eta_star1*N*0.5)/eta_star1-l);
            z2 = abs((X2+eta_star2*N*0.5)/eta_star2-l);
            [xx1,p1]=min(abs((X1+eta_star1*N*0.5)/eta_star1-l));
            [xx2,p2]=min(abs((X2+eta_star2*N*0.5)/eta_star2-l));
        end
        
        function value = Price(S1,S2,K,r,sigma1,sigma2,rho,tau)
            [p1,p2,C,H] = Spread.FFTMatrix(S1,S2,K,r,sigma1,sigma2,rho,tau);
            Vmat=K*exp(-0.05*tau)*real(C.*ifft2(H));
            value = Vmat(p1,p2);
        end
        
        function [delta1,delta2,theta,gamma11,gamma12,gamma22] = Greek(S1,S2,K,r,sigma1,sigma2,rho,tau)
            tic
            [p1,p2,C,HH1,HH2,HT,HHH11,HHH12,HHH22] = Spread.FFTMatrix(S1,S2,K,r,sigma1,sigma2,rho,tau);
            
            delta1mat=K/S1*exp(-0.05*tau)*real(C.*ifft2(HH1));
            delta2mat=K/S2*exp(-0.05*tau)*real(C.*ifft2(HH2));
            thetamat=K*exp(-0.05*tau)*real(C.*ifft2(HT));
            gamma11mat=-1/S1*delta1mat+K/S1^2*exp(-0.05*tau)*real(C.*ifft2(HHH11));
            gamma12mat=K/(S1*S2)*exp(-0.05*tau)*real(C.*ifft2(HHH12));
            gamma22mat=-1/S2*delta2mat+K/S2^2*exp(-0.05*tau)*real(C.*ifft2(HHH22));
            p1
            p2
            delta1 = delta1mat(p1,p2);
            delta2 = delta2mat(p1,p2);
            theta=thetamat(p1,p2);
            gamma11 = gamma11mat(p1,p2);
            gamma12 = gamma12mat(p1,p2);
            gamma22 = gamma22mat(p1,p2);
            toc
        end
        
        function u_test1 = utest(u_bar,X0,N)           
            for i=0:N-1
                u_test1 = pi * (i - N/2) / X0;
                if u_test1 > u_bar
                    return
                end
            end
        end
        
        function value = P_hat(u)
            value = cgamma(1i*(u(1)+u(2))-1) * cgamma(-1i*u(2)) / cgamma(1i*u(1)+1);
        end
        
        function value = Phi(mu,Sigma,tau)
            value = exp(1i*mu*tau - Sigma*tau*0.5);
        end
        
        function Plot(K,r,sigma1,sigma2,rho,tau,M)
            N=M/25;
            x = N:N:M;
            X = meshgrid(x,x);
            graphdelta1=zeros(M/N);
            graphdelta2=zeros(M/N);
            graphtheta=zeros(M/N);
            graphgamma11=zeros(M/N);
            graphgamma12=zeros(M/N);
            graphgamma22=zeros(M/N);
            
            for j=1:M/N
                for i=1:M/N
                    S1=X(i,j);
                    S2=X(j,i);
                    [delta1,delta2,theta,gamma11,gamma12,gamma22] = Spread.Greek(S1,S2,K,r,sigma1,sigma2,rho,tau);
                    graphdelta1(i,j)=delta1;
                    graphdelta2(i,j)=delta2;
                    graphtheta(i,j)=theta;
                    graphgamma11(i,j)=gamma11;
                    graphgamma12(i,j)=gamma12;
                    graphgamma22(i,j)=gamma22;
                end
            end
            
            figure()
            surf(N:N:M,N:N:M,graphdelta1)
            zlim([0 1])
            xlabel('Price of Asset 1')
            ylabel('Price of Asset 2')
            zlabel('Delta 1')
            figure()
            surf(N:N:M,N:N:M,graphdelta2)
            zlim([-1 0])
            xlabel('Price of Asset 1')
            ylabel('Price of Asset 2')
            zlabel('Delta 2')
            figure()
            surf(N:N:M,N:N:M,graphtheta)
%             zlim([0 0.1])
            xlabel('Price of Asset 1')
            ylabel('Price of Asset 2')
            zlabel('Theta')
            figure()
            surf(N:N:M,N:N:M,graphgamma11)
            zlim([0 0.1])
            xlabel('Price of Asset 1')
            ylabel('Price of Asset 2')
            zlabel('Gamma11')
            figure()
            surf(N:N:M,N:N:M,graphgamma12)
            zlim([-0.1 0])
            xlabel('Price of Asset 1')
            ylabel('Price of Asset 2')
            zlabel('Gamma12')
            figure()
            surf(N:N:M,N:N:M,graphgamma22)
            zlim([0 0.1])
            xlabel('Price of Asset 1')
            ylabel('Price of Asset 2')
            zlabel('Gamma22')
        end
    end
end

