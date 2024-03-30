function [X,Y] = loaddata(ind)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
if ind == 1
    load('glass_uni.mat');
elseif ind ==2
    load('yeast_uni.mat');
elseif ind ==3
    load('wine_uni.mat');
elseif ind ==4
    load('german_uni.mat');
elseif ind ==5
    load('dermatology_uni.mat');
elseif ind ==6
    load('JAFFE_1024.mat');
end
X = double(X);
Y = Y;
end


