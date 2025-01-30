function [m]=tomas(A,d)
%A is sparse
%solve linear equation system, A is a matrix and d is vector of results and
%m is the answer vector
I=length(d); %get the length of the vector of results
m=zeros(1,I);%preallocation
%fix the matrix for 3- diagonal matrix
d(I)=d(I)-d(I-2)*(A(I,I-3)/(A(I-2,I-3)));
A(I,:)=A(I,:)-A(I-2,:).*(A(I,I-3)/(A(I-2,I-3)));
d(I)=d(I)-d(I-1)*(A(I,I-2)/(A(I-1,I-2)));
A(I,:)=A(I,:)-A(I-1,:).*(A(I,I-2)/(A(I-1,I-2)));
A(I,:)=A(I,:)-A(I-2,:).*(A(I,I-3)/(A(I-2,I-3)));
A(I,:)=A(I,:)-A(I-1,:).*(A(I,I-2)/(A(I-1,I-2)));

E=zeros(1,I);%preallocation
F=E;%preallocation
E(1)=A(1,2)/A(1,1); F(1)=d(1)/A(1,1);%intial condition
for L=2:I %run the algorithm
    if(L~=I)
       E(L)=(A(L,L+1))/(A(L,L)-E(L-1)*A(L,L-1)); 
    end
    F(L)=(d(L)-F(L-1)*A(L,L-1))/(A(L,L)-E(L-1)*A(L,L-1));
end
%applied the solutions:
m(I)=F(I);
for l=I-1:-1:1
    m(l)=F(l)-E(l)*m(l+1);
end
end