
%%% Matlab program for error correction scheme for fully correlated channels on n qubits
%
%%% Speicify an integer n > 1, the following commands will generate the encoding matrix P_n
%
%%% CK Li May 2019
if mod(n,2) == 1   
       P = eye(8); P3 = P(:,[1,6,4,7,8,3,5,2]); Pn = P3;
          k = (n-1)/2;
       for j = 2:k
            Pn = kron(eye(4),Pn)*kron(P3,eye(2^(2*j-2)));
   end 
else 
     H = [1 1; 1 -1]/sqrt(2); C01 = [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]; 
                    P2 = C01*kron(eye(2),H)*C01;   
   if n == 2
              Pn = P2;
   else
       P = eye(8); P3 = P(:,[1,6,4,7,8,3,5,2]); Pn = P3; k = (n-2)/2;
     for j = 2:k
           Pn = kron(eye(4),Pn)*kron(P3,eye(2^(2*j-2)));
     end  
           Pn = kron(eye(2),Pn)*kron(P2,eye(2^n/4));
   end
end
%
%   Set up the error operators for the channel
%
     X = [0 1; 1 0];  Y = [0 -i; i 0];   Z = [1 0; 0 -1];  
         Xn = X; Yn = Y; Zn = Z;
   for j = 2:n
        Xn = kron(X,Xn); Yn = kron(Y,Yn); Zn = kron(Z,Zn);
   end
% 
%%% Test the encoding and decoding schemes using the following commands.
%
%%% If n is odd, check
%
%     (Pn'*Xn*Pn, Pn'*Yn*Pn, Pn'*Zn*Pn) = (kron(X,I), (-1)^k*Kron(Y,I), Kron(Z,I))
%
    II = eye(2^(n-1));  norm(Pn'*Xn*Pn -kron(X,II)), 
    norm(Pn'*Yn*Pn - (-1)^k*kron(Y,II)), norm(Pn'*Zn*Pn -kron(Z,II))
%
%  Next, we verify 
%
%       Pn'(En(Pn(A\otimes B)Pn')Pn = \hat A \otimes B for A in D_2, B in D_{2^n/2}
%
%  Generate random S in D_2, S in D_{2k}
%
   S = rand(2,2) + i*rand(2,2) - rand(1,1)*(1+i)*eye(2); S = S*S'; S = S/trace(S); 
   K = 2^(n-1); 
   R = rand(K,K) + i*rand(K,K) - rand(1,1)*(1+i)*eye(K); R = R*R'; R = R/trace(R);
%
%  Encode kron(S,R) and compared with the decoded state for each error operator.
%
   A = kron(S,R); AA = Pn*A*Pn'; 
   norm( Pn'*Xn*AA*Xn*Pn - kron(X*S*X',R)), norm(Pn'*Yn*AA*Yn*Pn - kron(Y*S*Y',R))
   norm( Pn'*Zn*AA*Zn*Pn - kron(Z*S*Z',R))
%
%%% If n is even, check
%
%        (Pn'*Xn*Pn, Pn'*Yn*Pn, Pn'*Zn*Pn) = (kron(Dx,I), Kron(Dy,I), Kron(Dz,I))
%
       Dx = diag([1 -1 1 -1]); Dy = diag([-1 -1 1 1]); Dz = diag([1 -1 -1 1]);
       II = eye(2^n/4);  norm(Pn'*Xn*Pn - kron(Dx,II)); 
     norm(Pn'*Yn*Pn - (-1)^k*kron(Dy,II)), norm(Pn'*Zn*Pn - kron(Dz,II))
%
%%% Next, we verify if A corresponds to 2 classical bits, and B is (n-2) qubits.
%
%     Pn'(En(Pn(A\otimes B)Pn')Pn = A \otimes B for A in D_2, B is D_{2^n/2}
%

%  Set up the classical bits in D_4, and arbitrary qubits in D_{2k}
%
   K = 2^(n-2); 
   R = rand(K,K) + i*rand(K,K) - rand(1,1)*(1+i)*eye(K); R = R*R'; R = R/trace(R);
   b0 = [1 0; 0 0]; b1 = [0 0; 0 1]; 
   b00 = kron(b0,b0); b01 = kron(b0,b1); b10 = kron(b1,b0); b11 = kron(b1,b1);
%
%  Encode and decode
%
   S = b00; A = kron(S,R); AA = Pn*A*Pn'; 
   norm(Pn'*Xn*AA*Xn*Pn - kron(S,R)), norm(Pn'*Yn*AA*Yn*Pn - kron(S,R))
   norm(Pn'*Zn*AA*Zn*Pn - kron(S,R))
%
   S = b01; A = kron(S,R); AA = Pn*A*Pn'; 
   norm(Pn'*Xn*AA*Xn*Pn - kron(S,R)), norm(Pn'*Yn*AA*Yn*Pn - kron(S,R))
   norm(Pn'*Zn*AA*Zn*Pn - kron(S,R))
%
   S = b10; A = kron(S,R); AA = Pn*A*Pn'; 
   norm(Pn'*Xn*AA*Xn*Pn - kron(S,R)), norm(Pn'*Yn*AA*Yn*Pn - kron(S,R))
   norm(Pn'*Zn*AA*Zn*Pn - kron(S,R))
%
   S = b11; A = kron(S,R); AA = Pn*A*Pn'; 
   norm(Pn'*Xn*AA*Xn*Pn - kron(S,R)), norm(Pn'*Yn*AA*Yn*Pn - kron(S,R))
   norm(Pn'*Zn*AA*Zn*Pn - kron(S,R))
%%%
