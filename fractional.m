%% 返回分数t时刻的y(:,:) 
% y_pre是分数
function [y]=fractional(y_pre,t);
    global I J a d lambda sigma C_t b c;
    cvx_begin quiet
        cvx_solver Mosek;
        variable x(J) nonnegative;
        variable y(I,J) nonnegative;
        variable z(J) nonnegative;
        variable w nonnegative;

        %% 定义问题

        p_t1=sum(a * sigma(:,t) .* x,1);%一维

        temp2 = b(:,t) * sigma(:,t)' .* y;
        p_t2=sum(temp2(:)); %二维

        epsilon=0.1;
        delta=log(1+(1/epsilon));
        %pre是y的前一个t的结果 和y同形

        %temp3 = c / delta .* ( rel_entr(y+epsilon,y_pre+epsilon) + y_pre - y);
%         temp3 = c / delta .* ( rel_entr(y+epsilon,y_pre+epsilon) - y + y_pre);
        temp3 = c / delta .* ( rel_entr(y+epsilon,y_pre+epsilon) - y);
        p_t3 = sum(temp3(:)); %二维

        p_t4=d*w;

        p_t = p_t1 + p_t2 + p_t3 + p_t4;


        %% 解

        minimize(p_t)
        subject to
            %1a
            x(:)+sum(y(:,:),1)' >= lambda(:,t);
            %1b
            x()+sum(y(:,:),1)' <= 1;
            %1c
            x>=z;
            %1d(2b)
            C_t(t) * z + sigma(:,t) .* sum(y(:,:),1)' == sigma(:,t);
            %1e 
            w + sum(z(:),1) == 1;
    cvx_end
end
