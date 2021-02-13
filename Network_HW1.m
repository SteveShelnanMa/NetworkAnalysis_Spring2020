
%%   5)
% A)
% Construct the binary adjacency matrix A for digraph
A = zeros(size(A,1),size(A,2))

for i = 1:size(A,1)
    for j = 1:size(A,2)
        A(i,j) = (sum(Y(i,j,:))>0)
    end
end

% B)
% Compute the number of directed edges (arcs) in the network
sum(A,"all")
issymmetric(A)

% C) 
% firstly construct adjacency matrix B for undirectred graph
B = zeros(size(A,1),size(A,2))

for i = 1:size(A,1)
    for j = 1:size(A,2)
        B(i,j) = (A(i,j)+A(j,i)>0)
    end
end

issymmetric(B)
% then compute the number of undirected edges in the network
sum(B,"all")/2


% D)
% Compute the number of mutual arcs in the network
C = A .* transpose(A)

sum(C,"all")/2


% E) 
n = 0
for i = 1:size(A,1)
    if sum(A(:,i)) == 0
        n = n + 1
        disp(employees(i))
    end
end


% F)
n = 0
for i = 1:size(A,1)
    if sum(A(i,:)) == 0
        %n = n + 1
        disp(employees(i))
    end
end


% G)
n = 0
for i = 1:size(A,1)
    if sum(A(:,i)) >= 30
        n = n + 1    
    end
end


% H)
n = 0
for i = 1:size(A,1)
    if sum(A(i,:)) >= 30
        n = n + 1    
    end
end

% I)
isdirected = 1

triangle(A,isdirected)

% J)
% construct the degree matrix D
D = zeros(size(B,1),size(B,2))

for i = 1:size(B,1)
    D(i,i) = sum(B(i,:))
end

issymmetric(D)

% construct the Laplacian of the undirected graph
L = D - B
% get the eigenvalue of L
e = eig(L)
% get the number of zero eigenvalues of L
sum(e<10^-13)

%% 6)
% A)
A_Fig1 = [0 1 1 0 0 0 0;
          1 0 1 1 0 0 0;
          1 1 0 0 1 0 0;
          0 1 0 0 1 1 1;
          0 0 1 1 0 1 0;
          0 0 0 1 1 0 1;
          0 0 0 1 0 1 0]

issymmetric(A_Fig1)  
      
D_Fig1 = zeros(size(A_Fig1,1),size(A_Fig1,2))

for i = 1:size(A_Fig1,1)
    D_Fig1(i,i) = sum(A_Fig1(i,:))
end     
      
L_Fig1 = D_Fig1 - A_Fig1     
 
% check the two definitions are equivalent
L1 = D_Fig1^(-.5)*L_Fig1*D_Fig1^(-.5)     
L2 = eye(size(A_Fig1,1)) - D_Fig1^(-.5)*A_Fig1*D_Fig1^(-.5) 
L1 - L2 

% check the Laplacian and normalized Laplacian have non_negative
% eigenvalues
eig(L_Fig1)
eig(L1)

% check the eigenvectors
L_Fig1*([1,1,1,1,1,1,1].')

L1*D_Fig1^(0.5)*([1,1,1,1,1,1,1].')
L2*D_Fig1^(0.5)*([1,1,1,1,1,1,1].')

% B)

[eigvec,eigval] = eig(L_Fig1)
[eigvec,eigval] = eig(L1)

% C) modify the Fig 1

A_Fig1 = [0 1 1 0 0 0 0;
          1 0 1 1 0 0 0;
          1 1 0 1 0 0 0;
          0 1 1 0 0 0 0;
          0 0 0 0 0 1 1;
          0 0 0 0 1 0 0;
          0 0 0 0 1 0 0]


issymmetric(A_Fig1)  
      
D_Fig1 = zeros(size(A_Fig1,1),size(A_Fig1,2))

for i = 1:size(A_Fig1,1)
    D_Fig1(i,i) = sum(A_Fig1(i,:))
end     
      
L_Fig1 = D_Fig1 - A_Fig1     
 
% check the two definitions are equivalent
L1 = D_Fig1^(-.5)*L_Fig1*D_Fig1^(-.5)     
L2 = eye(size(A_Fig1,1)) - D_Fig1^(-.5)*A_Fig1*D_Fig1^(-.5) 
L1 - L2 

% check the Laplacian and normalized Laplacian have non_negative
% eigenvalues
eig(L_Fig1)
eig(L1)

% check the eigenvectors
L_Fig1*([1,1,1,1,1,1,1].')

L1*D_Fig1^(0.5)*([1,1,1,1,1,1,1].')
L2*D_Fig1^(0.5)*([1,1,1,1,1,1,1].')


[eigvec,eigval] = eig(L_Fig1)
[eigvec,eigval] = eig(L1)

 







