function opts = loadconfig(ind)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
if ind == 1
    opts.anchor_rate = 0.2;
    opts.order = 4;
    opts.anchorSelectStyle = 3;
    rng(8);
elseif ind ==2
    opts.anchor_rate = 0.1;
    opts.order = 2;
    opts.anchorSelectStyle = 1;
elseif ind ==3
    opts.anchor_rate = 0.1;
    opts.order = 3;
    opts.anchorSelectStyle = 3;
    rng(10);
elseif ind ==4
    opts.anchor_rate = 0.1;
    opts.order = 2;
    opts.anchorSelectStyle = 3;
    rng(10);
elseif ind ==5
    opts.anchor_rate = 0.1;
    opts.order = 2;
    opts.anchorSelectStyle = 3;
    rng(3);
elseif ind ==6
    opts.anchor_rate = 0.3;
    opts.order = 5;
    opts.anchorSelectStyle = 3;
    rng(1);
end
end

