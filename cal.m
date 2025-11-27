function [f,f1,f2,f3,f4]=cal(x,y,w);
    global I J T sigma a b d c;
    f = zeros(T,1);
    f1 = zeros(T,1);
    f2 = zeros(T,1);
    f3 = zeros(T,1);
    f4 = zeros(T,1);
    for t=1:T
        %p_t = p_t1 + p_t2 + p_t3 + p_t4;
        p_t1=sum(a .* sigma(:,t) .* x(:,t),1);%一维

        temp2 = b(:,t) * sigma(:,t)' .* y(:,:,t);
        p_t2=sum(temp2(:)); %二维

        if t==1
            yy=y(:,:,t);
        else
            yy=y(:,:,t)- y(:,:,t-1); %处理()^+
        end
        for i=1:I
            for j=1:J
                if yy(i,j)<0 
                   yy(i,j)=0;
                end
            end
        end
        temp3=c .* yy;
        p_t3 = sum(temp3(:));

        p_t4=d*w(t);

        f(t,1) = p_t1 + p_t2 + p_t3 + p_t4;
        f1(t,1) = p_t1;
        f2(t,1) = p_t2;
        f3(t,1) = p_t3;
        f4(t,1) = p_t4;
    end
end