%% y和y_pre都是整数(:,:,t)
% x是分数
function [z,w]=fractional3(x,y,y_pre,t);
    global I J a d lambda sigma C_t b c;
    cvx_begin quiet
        cvx_solver Mosek;
        variable z(J) nonnegative;
        variable w nonnegative;

        %% 定义问题

        p_t1=sum(a * sigma(:,t) .* x,1);%一维

        temp2 = b(:,t) .* sigma(:,t)' .* y;
        p_t2=sum(temp2(:)); %二维

        epsilon=0.1;
        delta=log(1+(1/epsilon));
        %pre是y的前一个t的结果 和y同形

        %temp3 = c / delta .* ( rel_entr(y+epsilon,y_pre+epsilon) + y_pre - y);
%         temp3 = c / delta .* ( rel_entr(y+epsilon,y_pre+epsilon) - y + y_pre);
        temp3 = c / delta .* ( rel_entr(y+epsilon,y_pre+epsilon) - y );
        p_t3 = sum(temp3(:)); %二维

        p_t4=d*w;

        p_t = p_t1 + p_t2 + p_t3 + p_t4;


        %% 解

        minimize(p_t)
        subject to
            %1a
            x+sum(y,1)' >= lambda(:,t);
            %1b
            x+sum(y,1)' <= 1;
            %1c
            x>=z;
            %1d(2b)
            C_t(t) * z + sigma(:,t) .* sum(y,1)' == sigma(:,t);
            %1e 
            w + sum(z,1) == 1;
%         for j=1:J
%             %1a
%             x(j)+sum(y,1) >= lambda(j,t);
%             %改 1b
%             x(j)+sum(y,1) <= 1;
%             %1c
%             x(j)>=z(j);
%             %1d(2b)
%             C_t(t) * z(j) + sigma(j,t) * sum(y(:,j),1) >= sigma(j,t);
%             %1e (2c, M)
%             w + sum(z) == 1;
%             %2a 
%             %sum(x)-x(j)+sum(y(:))-sum(y,1)>=sum(lambda(:,t))-1
%         end

    cvx_end
end
