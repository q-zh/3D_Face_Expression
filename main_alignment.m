
clear;
addpath(genpath('depend'));
data_folder = './data1/';
process_data;

Options.normalWeighting = 0;
Options.plot = 1;

our_output = in_template;
tic
[pointsTransformed, X] = our_nricp(in_template, input, Options, trans_para.index1(1:end), trans_para.index2(1:end));
toc
[our_output.vertices, landmark] = transform_back(pointsTransformed, trans_para);
figure;ShowFace(our_output);view(180,90);

nricp_output = in_template;
tic
[pointsTransformed, X] = nricp(in_template, input, Options, trans_para.index1(1:end), trans_para.index2(1:end));
toc
[nricp_output.vertices, landmark] = transform_back(pointsTransformed, trans_para);
figure;ShowFace(nricp_output);view(180,90);