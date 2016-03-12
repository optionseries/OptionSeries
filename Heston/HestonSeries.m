function [vsum, vmat]=HestonSeries(u,S,K,t,r,theta,v,rho,kappa,eta,M)
x=log(S/K);
t=t./M;
r=r.*M;
v=v.*M;
theta=theta.*M;
kappa=kappa.*M;
eta=eta.*M;

a=1/2-r/theta;
b=zeros(8,8);
for i=1:8
    for j=1:8
        d0=cell2mat(u(i,j));
        c1=d0(:,1);
        c2=a.^d0(:,2);
        c3=t.^d0(:,3);
        c4=x.^d0(:,4);
        c5=theta.^d0(:,5);
        c6=rho.^d0(:,6);
        c7=kappa.^d0(:,7);
        c8=exp(d0(:,8).*kappa.*t);
        
        b(i,j)=sum(c1.*c2.*c3.*c4.*c5.*c6.*c7.*c8).*(eta/(eta+1))^(i-1).*((v-theta)/(v-theta+1))^(j-1);
    end
end

b=b.*exp(a.*x-theta*(1-a)^2*t/2-x.^2/(2.*theta.*t))/sqrt(8.*pi.*theta.*t);
b(1,1)=exp(x).*erfc(-(x+(1-a).*theta.*t)/sqrt(2.*theta.*t))/2-exp(-r.*t).*erfc(-(x-a.*theta.*t)/sqrt(2.*theta.*t))/2;
b=b.*K;
vmat=b;
vsum=sum(b(:));
end