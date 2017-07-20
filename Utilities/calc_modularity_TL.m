function [Q Qds Qds_Li C] = calc_modularity_TL(Ci,B)
% Calculate different quality metrics

 K=sum(B);                               %degree
 N=length(B);                            %number of vertices
 m=sum(K);  
% B=B-(K.'*K)/m;  
% 
% s=Ci(:,ones(1,N)); 
% 
% Q=~(s-s.').*B/m;
% Q=sum(Q(:));    
%     

%

G = unique(Ci);
G(G==-1) = []; % Ignore nodes assigned alone

for i = 1:length(G)
      gin = Ci==G(i);
      gout = Ci~=G(i);
      ginnum = nnz(gin);
      Ein = sum(sum(B(logical(gin),logical(gin))));
      Eout = sum(sum(B(logical(gin),logical(gout))));
      
      Gmod(i) = Ein/m - ((Ein+Eout)/m).^2;
      
      dci = Ein/(ginnum*ginnum);
      
      for j =1:length(G)
          
          gin2 = Ci==G(j);
          gin2num = nnz(gin2);
          Ecicj = sum(sum(B(logical(gin),logical(gin2))));
          dcicj = Ecicj./(ginnum*gin2num);
          dense_sub(j) = (Ecicj/m).*dcicj;
      end
      
      dense_sub_sum = sum(dense_sub);
      
      Gmod_dense(i) = Ein/m - ((Ein+Eout)/m).^2 - dense_sub_sum; % Chen, et al 2012, A New Metric for Quality of Network Community Structure
      
      Gmod_dense_Li(i) = (Ein-Eout)/ginnum; % From Li, et al 2008, Quantitative function for community detection
      
      C(i) = Eout/(Ein+Eout);
end

Q = sum(Gmod);
Qds = sum(Gmod_dense);
Qds_Li = sum(Gmod_dense_Li);
C = sum(C);
