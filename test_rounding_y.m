function [y]=test_rounding_y(y)
    global I J;
    for j=1:J
        inti=[];
        doublei=[];
        for i=1:I
            if (y(i,j)==0 || y(i,j)==1)
                inti=[inti i]; %只用记i就行 整数下标
            else
                doublei=[doublei i]; %分数下标
            end
        end
%         before=sum(y,1);
        while length(doublei)>1
            Iran=(randperm(length(doublei))); %randperm随机打乱一个序列 % 1～length的序列
            i1=doublei(Iran(1));
            i2=doublei(Iran(2));
%             if Iran(2)<Iran(1)
%                 temp=i2;
%                 i2=i1;
%                 i1=i2;
%             end
            rho1=min(1-y(i1,j),y(i2,j));
            rho2=min(y(i1,j),1-y(i2,j));
            if rand<(rho2/(rho1+rho2))
                y(i1,j)=y(i1,j)+rho1;
                y(i2,j)=y(i2,j)-rho1;
            else
                y(i1,j)=y(i1,j)-rho2;
                y(i2,j)=y(i2,j)+rho2;            
            end
            if (y(i2,j)==0||y(i2,j)==1)
                doublei(Iran(2))=[];
%             end
            elseif (y(i1,j)==0||y(i1,j)==1)
                doublei(Iran(1))=[];
            end        
        end
        if length(doublei)==1
            i1=doublei(1);
            if y(i1,j) < 0.01
                y(i1,j)=0;
            else 
                y(i1,j)=1;
            end
        end
    end
end