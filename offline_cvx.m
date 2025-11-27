%% 返回分数t时刻的y(:,:) 
% y_pre是分数
function [f,x,y,z,w,u]=offline_cvx();
    global I J T a d lambda sigma C_t b c;
    cvx_begin quiet
        cvx_solver Gurobi_2;
        variable x(J,T) binary;
        variable y(I,J,T) binary;
        variable u(I,J,T) nonnegative;
        variable z(J,T) nonnegative;
        variable w(T,1) nonnegative;

        %% 定义问题

        temp1=a * sigma .* x;
        p1=sum(temp1(:));
        
        temp2 = b * sigma';
        temp2=repmat(temp2,1,1,T);
        temp2=temp2.* y;
        p2=sum(temp2(:));
     
        temp3=repmat(c,1,1,T);
        temp3=temp3.*u;
        p3 = sum(temp3(:));

        temp4=d*w;
        p4=sum(temp4(:));
        
        p = p1 + p2 + p3 + p4;

        %% 解

        minimize(p)
        subject to
        %1a
        x+reshape(sum(y(:,:,:),1),J,T) >= lambda;
        %1b
        x+reshape(sum(y(:,:,:),1),J,T) <= 1;
        %1c
        x>=z;
        %1d(2b)
        C_t_r=repmat(C_t,J,1);    
        C_t_r .* z + sigma .* reshape(sum(y,1),J,T) == sigma;
        %1e 
        w' + sum(z,1) == 1;      
        u(:,:,1)>=y(:,:,1);
        u(:,:,2:T)>=y(:,:,2:T)-y(:,:,1:T-1);           
        
%         for j=1:J
%             for t=1:T
%                 %1a
%                 x(j,t)+sum(y(:,j,t),1) >= lambda(j,t);
%                 %1b
%                 x(j,t)+sum(y,1) <= 1;
%                 %1c
%                 x(j,t)>=z(j,t);
%                 %1d(2b)
%                 C_t(t) * z(j,t) + sigma(j,t) * sum(y(:,j,t),1) >= sigma(j,t);
%                 %1e 
%                 w(t) + sum(z(:,t),1) == 1;
%                 for i=1:I
%                     if t==1
%                         u(i,j,t)>=y(i,j,t);
%                     else
%                         u(i,j,t)>=y(i,j,t)-y(i,j,t-1);
%                     end
%                 end
%             end
%         end
    cvx_end
    p=cvx_optval;
    [f,f1,f2,f3,f4]=cal(x,y,w); 
end
