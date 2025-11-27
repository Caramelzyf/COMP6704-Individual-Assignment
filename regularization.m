%overview
obj_f_test = [];
%正则化+round 
obj_f = [];
obj_x = zeros(J,T);
obj_y = zeros(I,J,T);
obj_z = zeros(J,T);
obj_w = zeros(T,1);

%正则化 不round
obj_frac_x = zeros(J,T);
obj_frac_y = zeros(I,J,T);

fractional_f=zeros(T);
for t=1:T
    %用t-1的分数y 算t的分数y
    if t==1
        obj_frac_y(:,:,t) = fractional(zeros(I,J),t);
    else
        obj_frac_y(:,:,t) = fractional(obj_y(:,:,t-1),t);
    end
    
    %round y
    obj_y(:,:,t)=test_rounding_y(obj_frac_y(:,:,t));
    
%     %用t-1的整数y和t的整数y 代回fractional2 求分数x
%     if t == 1
%         obj_frac_x(:,t)=fractional2(obj_y(:,:,t),zeros(I,J),t);
%     else
%         obj_frac_x(:,t)=fractional2(obj_y(:,:,t),obj_y(:,:,t-1),t);
%     end
    
    %用t的整数y round x
    if sum(obj_y(:,:,t),1)==1
        obj_x(:,t) = 0;
    else
        obj_x(:,t) = lambda(:,t);
    end
    
    %用t的整数x,整数y t-1的整数y 代回fractional3 求分数z w
    if t == 1
        [obj_z(:,t),obj_w(t)]=fractional3(obj_x(:,t),obj_y(:,:,t),zeros(I,J),t);
    else
        [obj_z(:,t),obj_w(t)]=fractional3(obj_x(:,t),obj_y(:,:,t),obj_y(:,:,t-1),t);
    end   
end
%算function的值 用obj xyzw
[obj_f,obj_f1,obj_f2,obj_f3,obj_f4]=cal(obj_x,obj_y,obj_w); 
