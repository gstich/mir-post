clear all;
clc;


path(path,'../');
% Run the sum planes scripts on the data and save the mean field.
ipath = '/p/lscratchrza/olson45/TBL/post_proc';
post = 'post';


% Coarse
%clear data;
t1 = 700;
tf = 900;
data = sum_planes([ipath,'med_nomu'],post,t1,tf);
save('../../data/jcp/med_nomu.mat','data');

t1 = 260;
tf = 460;
data = sum_planes([ipath,'med_mu'],post,t1,tf);
save('../../data/jcp/med_mu.mat','data');

t1 = 1000;
tf = 1196;
data = sum_planes([ipath,'med_muNew'],post,t1,tf);
save('../../data/jcp/med_muNew.mat','data');




% Medium
clear data;
t1 = 575;
tf = 775;
data = sum_planes([ipath,'fin_nomu'],post,t1,tf);
save('../../data/jcp/fin_nomu.mat','data');

t1 = 200;
tf = 350;
data = sum_planes([ipath,'fin_mu'],post,t1,tf);
save('../../data/jcp/fin_mu.mat','data');

t1 = 1100;
tf = 1400;
data = sum_planes([ipath,'fin_muNew'],post,t1,tf);
save('../../data/jcp/fin_muNew.mat','data');
