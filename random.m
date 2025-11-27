%纯随机
function [random_f,random_x,random_y,random_z,random_w]=random()
    global I J T a d lambda sigma C_t b c;
    random_x = zeros(J,T);
    random_y = zeros(I,J,T);
    random_z = zeros(J,T);
    random_w = zeros(T,1);
    capacity=C_t;
    for t=1:T
        for j=1:J
            if lambda(j,t)==1
                random=randi([0,I]); %生成0～I的随机数
                if random==0 && capacity(t)-sigma(j,t) >= 0
                    random_x(j,t)=1;
                    capacity(t)=capacity(t)-sigma(j,t);
                    random_y(:,j,t)=0;
                else
                    if random==0 %0但是本地容量不够
                        random=randi([1,I]);
                    end
                    random_y(random,j,t)=1;
                end    
            end
        end
        %固定x y 求zw
        cvx_begin quiet
            cvx_solver Mosek;
            variable z(J) nonnegative;
            variable w nonnegative;
            %% 定义问题
            p_t = d*w;
            %% 解
            minimize(p_t)
            subject to
                %1c
                random_x(:,t)>=z;
                %1d(2b)
                C_t(t) * z + sigma(:,t) .* sum(random_y(:,:,t),1)' == sigma(:,t);
                %1e 
                w + sum(z(:),1) == 1;
        cvx_end
        random_z(:,t) = z;
        random_w(t,1) = w;
    end    
    random_f=cal(random_x,random_y,random_w);
end