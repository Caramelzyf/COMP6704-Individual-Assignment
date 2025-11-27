%overview
obj_f_test = [];
%正则化+round 
%rfhc_x = zeros(J,T);
rfhc_y = zeros(I,J,T);
rfhc_z = zeros(J,T);
%rfhc_w = zeros(T,1);

%正则化 不round
%rfhc_frac_x = zeros(J,T);
rfhc_frac_y = zeros(I,J,T);

rfhc_f=zeros(T);
n=T/ww;
for nn=1:n %几轮
    
    %step1
    %一共ww个时间窗口 tend时保留正则化结果
    tstart=(nn-1)*ww+1;
    tend=nn*ww;    
    for t=tstart:tend
        %用t-1的分数y 算t的分数y
        if t==1
            rfhc_frac_y(:,:,t)= fractional_PFO(zeros(I,J),t);
        else
            rfhc_frac_y(:,:,t)= fractional_PFO(rfhc_frac_y(:,:,t-1),t);
        end         
%  rfhc_frac_y(:,:,t)= fractional(zeros(I,J),t);
    end   
    
    cvx_begin quiet
        cvx_solver Mosek;
        %variable x(J,ww-1) nonnegative;
        variable y(I,J,ww-1) nonnegative;
        variable u(I,J,ww) nonnegative;
        variable z(J,ww-1) nonnegative;
        %variable w(ww-1) nonnegative;

        %% 定义问题
        % p是sum t到t+w-1的问题且t+w-1时刻的值大部分都确定了
        sigma_r=sigma(:,tstart:tend-1);
        ct=C_t(1,tstart:tend-1);
        temp1 = a * sigma_r .* z;
        p_1=sum(temp1(:));

        b_r=reshape(b,size(b,1),1,size(b,2));
        b_r=repmat(b_r,1,J,1);
        b_r=b_r(:,:,tstart:tend-1);
        sigma_r=reshape(sigma,1,size(sigma,1),size(sigma,2));
        sigma_r=repmat(sigma_r,I,1,1);
        sigma_r=sigma_r(:,:,tstart:tend-1);
        temp2 = b_r .* sigma_r .* y;
        p_2=sum(temp2(:));

        c_r=reshape(c,size(c,1),size(c,2),1);
        c_r=repmat(c_r,1,1,ww);
        temp3 = c_r .* u;
        p_3 = sum(temp3(:));

        %temp4 = d * w;
        % p_4=sum(d*w,1);

        p = p_1 + p_2 + p_3; % + p_4;
        %% 解

        minimize(p)
        subject to
        %1a
        z+reshape(sum(y,1),J,[]) >= lambda(:,tstart:tstart+ww-2);
        %1b
        z+reshape(sum(y,1),J,[]) <= 1;
        %1c
        %x>=z; 
        %1d(2b)
        C_t_r=repmat(C_t,J,1);
        C_t_r(:,tstart:tstart+ww-2) .* z + sigma(:,tstart:tstart+ww-2) .* reshape(sum(y(:,:,:),1),J,[]) >= sigma(:,tstart:tstart+ww-2);
        %1e 
        %w' + sum(z,1) == 1;
        sum(z,1) <= 1;
        %u需要分成三段
        if tstart==1
            u(:,:,1)>=y(:,:,1);
        else
            u(:,:,1)>=y(:,:,1)-rfhc_frac_y(:,:,tstart-1);
        end
        u(:,:,2:ww-1)>=y(:,:,2:ww-1)-y(:,:,1:ww-2);
        u(:,:,ww)>=rfhc_frac_y(:,:,tend)-y(:,:,ww-1);
        
    cvx_end
    rfhc_frac_y(:,:,tstart:tend-1) = y; 

    for t=tstart:tend
        %round y
        rfhc_y(:,:,t)=test_rounding_y(rfhc_frac_y(:,:,t));


   
        % step2
        %用t的整数x,整数y t-1的整数y 代回fractional3 求分数z w
        if t == 1
            [rfhc_z(:,t)]=fractional3_PFO(rfhc_y(:,:,t),zeros(I,J),t);
        else
            [rfhc_z(:,t)]=fractional3_PFO(rfhc_y(:,:,t),rfhc_y(:,:,t-1),t);
        end 
    end
end
%算function的值 用obj xyzw
[rfhc_f,rfhc_f1,rfhc_f2,rfhc_f3]=cal_PFO(rfhc_y,rfhc_z);  

