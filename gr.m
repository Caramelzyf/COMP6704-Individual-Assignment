%%greedy+我们的round
%% 返回分数t时刻的y(:,:) 
% y_pre是分数
function [gr_f,gr_f1,gr_f2,gr_f3,gr_f4]=gr()
    global I J T a d lambda sigma C_t b c;
    gr_f=zeros(T,1);
    gr_x=zeros(J,T);
    gr_y=zeros(I,J,T);
    gr_w=zeros(T,1);
    
    %解xyzw的分数解
    for t=1:T
        cvx_begin quiet
            cvx_solver Mosek;
            variable x(J) nonnegative;
            variable y(I,J) nonnegative;
            variable u(I,J) nonnegative;
            variable z(J) nonnegative;
            variable w nonnegative;

            %% 定义问题
            temp1=a * sigma(:,t) .* x;%一维
            p1=sum(temp1(:));

            temp2 = b * sigma';
            temp2=temp2.* y;
            p2=sum(temp2(:)); %二维

            temp3=c.*u;
            p3 = sum(temp3(:)); %二维

            p4=d*w;

            p = p1 + p2 + p3 + p4;

            %% 解

            minimize(p)
            subject to
            %1a
            x(:)+sum(y(:,:),1)' >= lambda(:,t);
            %1b
            x()+sum(y(:,:),1)' <= 1;
            %1c
            x>=z;
            %1d(2b)
            C_t(t) * z + sigma(:,t) .* sum(y(:,:),1)' >= sigma(:,t);
            %1e 
            w + sum(z(:),1) == 1;
            if t==1
                u>=y;
            else
                u>=y-gr_y(:,:,t-1);
            end
        cvx_end
        gr_x(:,t)=x;
        gr_y(:,:,t)=y;
        
        %round y
        gr_y(:,:,t)=test_rounding_y2(gr_y(:,:,t));
        
        %求x
        if sum(gr_y(:,:,t),1)==1
            gr_x(:,t)=0;
        else
%             if rand<0.5
%                 gr_x(:,t)=lambda(:,t);
%             else
%                 gr_x(:,t)=1;
%             end
            gr_x(:,t)=lambda(:,t);
        end
        
        %固定xy 求zw
        cvx_begin quiet
            cvx_solver Mosek;
            variable z(J) nonnegative;
            variable w nonnegative;

            %% 定义问题
            p =d*w;

            %% 解
            minimize(p)
            subject to
            %1c
            x>=z;
            %1d(2b)
            C_t(t) * z + sigma(:,t) .* sum(gr_y(:,:,t),1)' == sigma(:,t);
            %1e 
            w + sum(z(:),1) == 1;
        cvx_end   
        gr_w(t,1)=w;        
    end
    [gr_f,gr_f1,gr_f2,gr_f3,gr_f4]=cal(gr_x,gr_y,gr_w); 
end
