clear
clc
addpath('F:\OXH\algorithm_fn_project\DeHaze')

f = imread('647_705-72-2.tif');
f= double(f);
% f=rot90(f,2);
record_max = max(f(:));
G = f./record_max;
a=5;
S0WinRow=a;
S0WinCol=a;
S0LT = CalculateInitialS0_LeftTop(G, S0WinRow, S0WinCol);% (G, WinRow, WinCol)
s_r = a-1; % s_r:= shift in row
s_c = a;% shift in col
ssc = round(s_c/2);
ssr = ssc-1;
S0LT1 = zeros(size(S0LT));
S0LT1(ssc:end,ssc:end) =  S0LT(1:end-ssr,1:end-ssr);
imwrite(uint16(65535*(S0LT1 )),['targetS0_S0LT1_Winsize',num2str(S0WinRow),'.tif']);
lambda=0.1;% default 0.01
StaticWinsize = 5;
[U, S, a_map, b_map] = TVGF_RemoveScattering(G, S0LT1, StaticWinsize, lambda);
imwrite(uint16(record_max*U),['targetU_S0LT1_Winsize',num2str(S0WinRow),'_StaticWinsize',num2str(StaticWinsize),'_lambda',num2str(lambda),'.tif']);
imwrite(uint16(record_max*S),['targetS_S0LT1_Winsize',num2str(S0WinRow),'_StaticWinsize',num2str(StaticWinsize),'_lambda',num2str(lambda),'.tif']);
imwrite(uint16(65535*a_map),['a_map_Winsize',num2str(S0WinRow),'_StaticWinsize',num2str(StaticWinsize),'_lambda',num2str(lambda),'.tif']);
imwrite(uint16(65535*b_map),['b_map_Winsize',num2str(S0WinRow),'_StaticWinsize',num2str(StaticWinsize),'_lambda',num2str(lambda),'.tif']);