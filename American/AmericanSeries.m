function p=AmericanSeries(d,u,S,K,t,r,v,N)
if v<r
    theta=r;
else
    theta=v;
end
A=-1/2-r/theta;
b=[];

for i=1:N
    d0=cell2mat(d(i));
    c1=d0(:,1);
    c2=A.^d0(:,2);
    c3=t.^d0(:,3);
    c4=(v-theta).^d0(:,4);
    c5=theta.^d0(:,5);
    
    b(end+1)=sum(c1.*c2.*c3.*c4.*c5);
end

if log(S/K)>sum(b)
    x=log(S/K)-sum(b);
else
    x=0;
end

v1=[];
v2=[];
for i=1:N
    u1=cell2mat(u(i,1));
    u2=cell2mat(u(i,2));
    c1=u1(:,1);
    c2=A.^u1(:,2);
    c3=t.^u1(:,3);
    c4=x.^u1(:,4);
    c5=(v-theta).^u1(:,5);
    c6=theta.^u1(:,6);
    
    c7=u2(:,1);
    c8=A.^u2(:,2);
    c9=t.^u2(:,3);
    c10=x.^u2(:,4);
    c11=(v-theta).^u2(:,5);
    c12=theta.^u2(:,6);
    
    v1(end+1)=sum(c1.*c2.*c3.*c4.*c5.*c6);
    v2(end+1)=sum(c7.*c8.*c9.*c10.*c11.*c12);
end

eab=exp((1/2-r/theta)*x-theta*A^2*t/2);
vexp=eab*exp(-x^2/(2*theta*t));
verfc=eab*erfc(x/sqrt(2*theta*t));

v=v1.*vexp+v2.*verfc;
if log(S/K)>sum(b)
    p=sum(v);
else
    p=sum(v)+exp(sum(b))-S/K;
end
p=K*p;
end