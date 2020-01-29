%% load data
[x,y,z,fl] = absload([data_folder '/data.abs']);
mask = (fl==1);
load([data_folder 'label.mat']);
load([data_folder 'landmark.mat']);
% figure; ShowFace(z,mask);

%% refine landmark
landmark = double(reshape(landmark(1,:,:), 68,2));
if landmark(49,1:2)==landmark(61,1:2)
    landmark(49,1) = landmark(49,1)-1;
end

if landmark(65,1:2)==landmark(55,1:2)
    landmark(55,1) = landmark(55,1)+1;
end

if landmark(62,1:2)==landmark(68,1:2)
    landmark(62,2) = landmark(62,2)-1;
    landmark(68,2) = landmark(68,2)+1;
end

if landmark(63,1:2)==landmark(67,1:2)
    landmark(63,2) = landmark(63,2)-1;
    landmark(67,2) = landmark(67,2)+1;
end

if landmark(64,1:2)==landmark(66,1:2)
    landmark(64,2) = landmark(64,2)-1;
    landmark(66,2) = landmark(66,2)+1;
end
landmark = landmark(18:end, :);

%% 2D landmark to 3D
mask_new = mask<-1;
for i = 1:size(landmark,1)
    mask_new(landmark(i,2), landmark(i,1)) = mask(landmark(i,2), landmark(i,1))==0;
end
[x,y] = meshgrid(1:640,1:480);
if sum(mask_new(:)) ~= 0
    z(mask_new) = griddata(x(mask),y(mask),z(mask),x(mask_new),y(mask_new));
    mask = or(mask, mask_new) ;
end
for i = 1:size(landmark,1)
    landmark(i,3)  = z(landmark(i,2), landmark(i,1));
end
% figure; ShowFace(z,mask, landmark);

%% crop face
nose_x = landmark(14,1);
nose_y = landmark(14,2);
nose_z = landmark(14,3);
dis_eye = sqrt(sum((abs((landmark(20,:)+landmark(23,:))/2 - (landmark(26,:)+landmark(29,:))/2)).^2));
l2_dis = sqrt((x-nose_x).^2 + (y-nose_y).^2 + (z-nose_z).^2);
mask_new = l2_dis < (dis_eye)*1.8; %dis_eye_h2*3.2; 1.2
mask = and(mask_new, mask);
% figure; ShowFace(z,mask, landmark);

%% crop mouth
mouth_points = landmark(32:43,:);
min_z = min(mouth_points(:,3));
box_left = min(mouth_points(:,1));
box_right = max(mouth_points(:,1));
box_up = min(mouth_points(:,2));
box_down = max(mouth_points(:,2));
z_mouth = z(box_up:box_down, box_left:box_right);
[g, d] = imgradient(z_mouth);
mask_mouth = and(g<50,z_mouth>min_z-10);
mask_moutht = mask==100;
mask_moutht(box_up:box_down, box_left:box_right) = ~mask_mouth;
in = inpolygon(x,y,double(landmark([44:end,44],1)),double(landmark([44:end,44],2)));
mask_mouth = or(mask_moutht, in);
% figure;ShowFace(vertex, edge);
%% load template
motion_array = ['an'; 'di'; 'fe'; 'ha'; 'sa'; 'su'; 'ne'];
t_a = [FER_result{1,1}.emotions.angry, FER_result{1,1}.emotions.disgust,...
    FER_result{1,1}.emotions.fear, FER_result{1,1}.emotions.happy, ...
    FER_result{1,1}.emotions.sad, FER_result{1,1}.emotions.surprise,...
    FER_result{1,1}.emotions.neutral];
label = motion_array(find(t_a==max(t_a(:))), :);
load(['template_' label '.mat']);
load('template_landmark.mat');
load('edgerlt.mat');
%% process data
in_template.faces = edgerlt';
in_template.vertices = template;
landmark_template = template(landmark_template,:);
landmark = transfer_landmark(landmark);
for i = 1:size(landmark_template, 1)
    temp  = sum((in_template.vertices - repmat(landmark_template(i,:), [size(in_template.vertices,1), 1])).^2,2);
    index1(i) = find(temp==min(temp(:)));
end
sign1 =  (landmark_template(1,2)+landmark_template(2,2))/2 - landmark_template(3, 2);
sign2 =  (landmark(1,2)+landmark(2,2))/2 - landmark(3, 2);
if sign1*sign2 < 0
    landmark(:,2) = 481 - landmark(:,2);
    landmark(:,1) = 641 - landmark(:,1);
    z = flip(z,1);
    mask = flip(mask,1);
    mask_mouth = flip(mask_mouth,1);
    z = flip(z,2);
    mask = flip(mask,2);
    mask_mouth = flip(mask_mouth,2);
end
[vertex, edge] = find_vertex_edge(z, mask, mask_mouth);
input.vertices = vertex;
input.faces = edge;
%% rough align
dis_eye1 = abs(landmark_template(1,1) - landmark_template(2,1));
dis_eye2 = abs(landmark(1,1) - landmark(2,1));
scale_eye = dis_eye1/dis_eye2;
input.vertices = input.vertices .* scale_eye;
landmark = landmark .* scale_eye;

% noise align
translate2 = landmark(3,:)-landmark_template(3,:);
input.vertices = input.vertices - repmat(translate2, [size(input.vertices,1),1]);
landmark = landmark - repmat(translate2, [size(landmark,1),1]);

% icp
R = [];
t = [];
[R, t] = icp(in_template.vertices(1:10:end,:)', input.vertices(1:500:end,:)');

input.vertices = R*input.vertices' + repmat(t, [1, size(input.vertices',2)]);
input.vertices = input.vertices';
landmark = R*landmark' + repmat(t, [1, size(landmark',2)]);
landmark = landmark';

% nose-xy align
translate1 = landmark(3,1:2)-landmark_template(3,1:2);
input.vertices(:, 1:2) = input.vertices(:,1:2) - repmat(translate1, [size(input.vertices,1),1]);
landmark(:, 1:2) = landmark(:, 1:2) - repmat(translate1, [size(landmark,1),1]);

% % scale up and dowm
dis_eye_h1 = abs((landmark_template(1,2) + landmark_template(2,2))/2 - landmark_template(3, 2));
dis_mouth_h1 = abs((landmark_template(5,2)+landmark_template(4,2))/2 - landmark_template(3, 2));

dis_eye_h2 = abs((landmark(1,2) + landmark(2,2))/2 - landmark(3, 2));
dis_mouth_h2 = abs((landmark(5,2)+landmark(4,2))/2 - landmark(3, 2));
scale_up =  dis_eye_h1/dis_eye_h2;
scale_down =  dis_mouth_h1/dis_mouth_h2;
nose_tip_y = landmark(3, 2);
input.vertices(find(input.vertices(:,2)>nose_tip_y), 2) = (input.vertices(find(input.vertices(:,2)>nose_tip_y), 2) - nose_tip_y).*scale_up + nose_tip_y;
input.vertices(find(input.vertices(:,2)<nose_tip_y), 2) = (input.vertices(find(input.vertices(:,2)<nose_tip_y), 2) - nose_tip_y).*scale_down + nose_tip_y;
landmark(find(landmark(:,2)>nose_tip_y), 2) = (landmark(find(landmark(:,2)>nose_tip_y), 2) - nose_tip_y).*scale_up + nose_tip_y;
landmark(find(landmark(:,2)<nose_tip_y), 2) = (landmark(find(landmark(:,2)<nose_tip_y), 2) - nose_tip_y).*scale_down + nose_tip_y;


for i = 1:size(landmark, 1)
    temp  = sum((input.vertices - repmat(landmark(i,:), [size(input.vertices,1), 1])).^2,2);
    index2(i) = find(temp==min(temp(:)));
end


trans_para.translate1 = translate1;
trans_para.translate2 = translate2;
trans_para.nose_tip_y = nose_tip_y;
trans_para.scale_eye = scale_eye;
trans_para.scale_up = scale_up;
trans_para.scale_down = scale_down;
trans_para.R = R;
trans_para.t = t;
trans_para.index1 = index1;
trans_para.index2 = index2;
trans_para.sign = sign1*sign2;

