clear;
addpath(genpath('depend'));
data_folder = './data2/';
process_data;

Options.normalWeighting = 0;
Options.plot = 1;
[pointsTransformed, X] = our_nricp(in_template, input, Options, trans_para.index1(1:end), trans_para.index2(1:end));

translate = mean(pointsTransformed) - mean(in_template.vertices);
in_rm = pointsTransformed - repmat(translate, [size(pointsTransformed,1),1]); 
result = expression_removal(in_template.vertices, label);
result =  result - in_template.vertices + pointsTransformed;

Output = in_template;
[Output.vertices, landmark] = transform_back(result, trans_para);
figure;ShowFace(Output);


nricp_output = in_template;
[nricp_output.vertices, landmark] = transform_back(pointsTransformed, trans_para);
figure;ShowFace(nricp_output);

img = imread([data_folder 'img.ppm']);
figure;imshow(img);
figure;img_fr = guide_2D_fr(img, nricp_output.vertices, Output.vertices, Output.faces);









