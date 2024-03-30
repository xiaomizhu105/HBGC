clc;clear all;
% 1-Glass; 2-Yeast; 3-Wine; 4-German; 5-Dermatology
% 6-JAFFE; 7-COIL-20; 8-MSRA-25
[X,Y] = loaddata(3); % 载入数据
opts = loadconfig(3);
% rng(10);

[n, dim] = size(X);
dimY = length(Y);
k = length(unique(Y));

if n ~=dimY
    X = X';
    dim = n;
    n = dimY;
end 

X = X./max(X,[],2);% n x dim
anchor_rate = opts.anchor_rate;% 锚点率
opts.IterMax = 50;
opts.anchorSelectStyle = opts.anchorSelectStyle; % 锚点选择方式
order = opts.order; % 高阶一共有几阶
gamma = 0.2;
m = fix(n*anchor_rate);
%% 锚点选择方式
disp('----------Anchor Selection----------');
if opts.anchorSelectStyle == 1 % direct sample
    [~,ind,~] = graphgen_anchor(X,m);
    centers = X(ind, :);% m x dim
elseif opts.anchorSelectStyle == 2 % rand sample
    vec = randperm(n);
    ind = vec(1:m);
    centers = X(ind, :);
elseif opts.anchorSelectStyle == 3 % KNP
    [~, ~, ~, ~, dis] = litekmeans(X,m);
    [~,ind] = min(dis,[],1);
    ind = sort(ind,'ascend');
    centers = X(ind, :);
elseif opts.anchorSelectStyle == 4 % kmeans sample
    [~, centers, ~, ~, ~] = litekmeans(X, m);
end
%% 一阶二部图构造
disp('----------1st order 2P Graphs Inilization----------');
D = L2_distance_1(X', centers'); % n x m
[~, idx] = sort(D, 2);
B1 = zeros(n,m);
for ii = 1:n
    id = idx(ii,1:k+1);
    di = D(ii, id);
    B1(ii,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
%% 高阶二部图构造
disp('----------Generate high order 2P Graphs----------');
B = cell(order,1); % cell;order x 1 -> n x m
BigLabmda = sum(B1,1);
B_n = B1.*(BigLabmda.^(-0.5));
% B_n = B1;

[U,sigma,Vt] = svd(B_n);

idx = randperm(n,m);
Um = U(idx,:);

for d = 1:order
    temp = U*(sigma*sigma').^d*Um';
    temp(temp<eps)=0;
    temp = temp./max(max(temp,[],2));
    %temp = temp./sum(temp,2);
    B{d,1} = temp;
end
% 
% %%
% disp('----------Optimization----------');
% alpha = 1/order*ones(order,1);% initial \alpha
% % initial P
% P = zeros(n,m);
% for o = 1:order
%     P = P + alpha(o)*B{o};
% end
% for iter = 1:opts.IterMax
%     fprintf('The %d-th iteration...\n',iter);
%     
%     % || Fix alpha and P,update F ||
%     [Fn,Fm,Dn,Dm] = max_Tr_FLsF(P,k);
%     
%     % || Fix alpha and F,update P ||
%     
%     % compute G
%     G = zeros(n,m);
%     for o = 1:order
%         G = G + 1/alpha(o)*B{o};
%     end
%     G = G/sum(1/alpha);
%     
%     %compute ti
%     fn = Dn*Fn;
%     fm = Dm*Fm;
%     dist = L2_distance_1(fn',fm');
%     
%     % update P
%     P = zeros(n,m);
%     for i=1:n
%         sub_term = G(i,:)-gamma/2*dist(i,:);
%         P(i,:) = EProjSimplex_new(sub_term);
%     end
%     
%     % compute sum(ev(1:k)) and sum(ev(1:k+1))
%     [sum_k,sum_k_1] = sum_ev(P,k);
%     if sum_k > 0.00000001
%         gamma = 2*gamma;
%     elseif sum_k_1 < 0.00000001
%         gamma = gamma/2;
%     else
%         fprintf('get k components P!!!\n');
%         break;
%     end  
%     
%     % || fix F,P, update \alpha ||
%     h = zeros(order,1);
%     for o = 1:order
%         h(o) = norm(P-B{o}, 'fro');
%     end
%     for o = 1:order
%         alpha(o) = h(o)/sum(h);
%     end
% end
% 
% ComputeAcc(P,Y);
% function ComputeAcc(P,Y)
%     [n, m] = size(P);
%     S=sparse(n+m,n+m);
%     S(1:n,n+1:end)=P; 
%     S(n+1:end,1:n)=P';
%     [c, y]=graphconncomp(S);
%     y1=y(1:n)';
%     i=1;
%     [result(i,:)] = ClusteringMeasure1(Y, y1);
%     fprintf('acc = %f !!\n',result(1));
% end







% B{1} = B1./max(max(B1,[],2));
% [U,sigma,Vt] = svd(B1);
% for d = 2:order
%     temp = U*sigma.^(2*d-1)*Vt';
%     temp(temp<eps)=0;
%     temp = temp./max(max(temp,[],2));
%     %temp = temp./sum(temp,2);
%     B{d,1} = temp;
% end
% 
%% 显示高阶二部图
[val,idx_per] = sort(Y(idx));
figure(1);
imagesc(B{1}(:,idx_per));
figure(2);
imagesc(B{2}(:,idx_per));
figure(3);
imagesc(B{3}(:,idx_per));

