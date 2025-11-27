global I;
global J;
global T;
global a;
global d;
global lambda;
global sigma;
global C_t;
global b;
global c;
global ww;


%% INTITIALIZA

T=100; % time slot
I=6; % scrubbing center
J=1000; % flow 
a=6.81; % operational cost
% d=0.84; % wasting cost
d=6500; % wasting cost

weight = 10;
lambda=load('D:/data/lambda.mat','-ascii');% lambda_jt 
sigma=load('D:/data/sigma.mat','-ascii');% sigma_jt 
% C_t=load('/Users/zy/Downloads/data/Ct.mat','-ascii');
C_t=load('D:/data/Ct.mat','-ascii');
b=load('D:/data/b.mat','-ascii');
c=load('D:/data/c_1000.mat','-ascii');

lambda=lambda(1:J,1:T); % lambda_jt 
sigma=sigma(1:J,1:T); % sigma_jt 
C_t=C_t(1,1:T);
% C_t=C_t .*3;
% r = 0.1 + 0.2 .* rand(1,T);
% C_t = C_t .* r;
% d=d.*100;
% C_t = C_t .* 3;

% d=d.*3;

b=b.*10; %为了平衡 b必须�?*10
% C_t=zeros(1,T);
b=b(1:I,1:T); % 
c=c(1:I,1:J); % 
% c=b * sigma' * weight;
c=c*weight;


% 引入noise
sum_b=sum(b,1);
sum_sigma=sum(sigma,1);
mean_b=sum_b./sum(lambda,1);
mean_sigma=sum_sigma./sum(lambda,1);

% % 准确的prediction 让p==0就行�?
noise_rfhc=[];%不准确的prediction 0 0.05 0.1 0.15
noise_opt=[];
ww=2; %预测窗口

ori_b=b;
ori_sigma=sigma;
[opt_f,opt_x,opt_y,opt_z,opt_w,opt_u]=offline_cvx();
noise_opt=[noise_opt,sum(opt_f)]
% for p = [0,0.05,0.1,0.15]
%     p
%     
%     noise_b=p*mean_b.*randn(I,T);
%     noise_sigma=p*mean_sigma.*randn(J,T);
%     b=ori_b+noise_b;
%     sigma=ori_sigma+noise_sigma;
%     lambda(sigma<=0)=0;
%     lambda(sigma>0)=1;
%     sigma(sigma<=0)=0;
%     
%     rfhc;
%     noise_rfhc=[noise_rfhc, sum(rfhc_f)]
% 
% end





% lambda=randi([0,1],[J,T]); % lambda_jt 
% sigma=randi([0,100],[J,T]); % sigma_jt 
% C_t=randi([5,2000],[1,T]);
% b=randi([1,10],[I,T]); % 
% weight=5;
% c=weight*b*sigma';

% tic
% regularization;
% sum_obj=sum(obj_f)
% toc
% % 

% tic
% rfhc;
% sum_rfhc=sum(rfhc_f)
% toc

tic
rfhc;
sum_rfhc=sum(rfhc_f)
toc
% 
% % 不管lambda(j,t)都设置xy
% tic
% [random_f,random_x,random_y,random_z,random_w]=random();
% sum_random=sum(random_f)
% toc
% 
% tic
% [gr_f,gr_f1,gr_f2,gr_f3,gr_f4]=gr();
% sum_gr=sum(gr_f)
% toc
% 
% tic
% [opt_f,opt_x,opt_y,opt_z,opt_w,opt_u]=offline_cvx();
% sum_opt=sum(opt_f)
% toc





