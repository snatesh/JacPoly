function K = promotion_mat(a,b,N)
NN = (0:(N-1))';
lam1 = (a+b+NN+1).*(a+b+NN+2)./((a+b+2*NN+1).*(a+b+2*NN+2)); lam1(1)=1;
lam2 = (a-b).*(a+b+NN+2)./((a+b+2*NN+2).*(a+b+2*NN+4));
lam3 = -(a+NN+2).*(b+NN+2)./((a+b+2*NN+4).*(a+b+2*NN+5));
K = spdiags([lam1,[0;lam2(1:N-1)],[0;0;lam3(1:N-2)]],0:2,N,N);
end
