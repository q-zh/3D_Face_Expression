function [test_face_neu_new, trans_landmark] = transform_back(test_face_neu, para)

translate1 = para.translate1;
translate2 = para.translate2; 
nose_tip_y = para.nose_tip_y; 
scale_eye = para.scale_eye;
scale_up = para.scale_up;
scale_down = para.scale_down;
R = para.R;
t = para.t;


trans_landmark = test_face_neu(para.index1, :);
trans_nose_tip_y = nose_tip_y;
test_face_neu_new = test_face_neu;
test_face_neu_new(find(test_face_neu_new(:,2)>trans_nose_tip_y), 2) = (test_face_neu_new(find(test_face_neu_new(:,2)>trans_nose_tip_y), 2) - nose_tip_y)./scale_up + trans_nose_tip_y;
test_face_neu_new(find(test_face_neu_new(:,2)<trans_nose_tip_y), 2) = (test_face_neu_new(find(test_face_neu_new(:,2)<trans_nose_tip_y), 2) - nose_tip_y)./scale_down + trans_nose_tip_y;
trans_landmark(find(trans_landmark(:,2)>trans_nose_tip_y), 2) = (trans_landmark(find(trans_landmark(:,2)>nose_tip_y), 2) - trans_nose_tip_y)./scale_up + trans_nose_tip_y;
trans_landmark(find(trans_landmark(:,2)<trans_nose_tip_y), 2) = (trans_landmark(find(trans_landmark(:,2)<nose_tip_y), 2) - trans_nose_tip_y)./scale_down + trans_nose_tip_y;

test_face_neu_new(:, 1:2) = test_face_neu_new(:,1:2) + repmat(translate1, [size(test_face_neu_new,1),1]);
trans_landmark(:, 1:2) = trans_landmark(:, 1:2) + repmat(translate1, [size(trans_landmark,1),1]);



% nVertsSource = size(test_face_neu_new, 1);
% I = (1:nVertsSource)';
% J = 4*I;
% D = sparse([I;I;I;I],[J-3;J-2;J-1;J],[test_face_neu_new(:);ones(nVertsSource,1)],nVertsSource, 4*nVertsSource);
% X = repmat([inv(R)'; -t'], nVertsSource, 1);
% test_face_neu_new =  D*X;

test_face_neu_new = inv(R) * (test_face_neu_new' - repmat(t, [1, size(test_face_neu_new,1)]));
trans_landmark = inv(R) * (trans_landmark' - repmat(t, [1, size(trans_landmark,1)]));
test_face_neu_new = test_face_neu_new';
trans_landmark = trans_landmark';
% nLandmark = size(landmark, 1);
% I = (1:nLandmark)';
% J = 4*I;
% D = sparse([I;I;I;I],[J-3;J-2;J-1;J],[landmark(:);ones(nLandmark,1)],nLandmark, 4*nLandmark);
% X = repmat([inv(R)'; -t'], nLandmark, 1);
% landmark =  D*X;

test_face_neu_new = test_face_neu_new + repmat(translate2, [size(test_face_neu_new,1),1]);
trans_landmark = trans_landmark + repmat(translate2, [size(trans_landmark,1),1]);

test_face_neu_new = test_face_neu_new ./ scale_eye;
trans_landmark = trans_landmark ./ scale_eye;


if para.sign<0
    test_face_neu_new(:,2) = 481 - test_face_neu_new(:,2);
    test_face_neu_new(:,1) = 641 - test_face_neu_new(:,1);
    
    trans_landmark(:,2) = 481 - trans_landmark(:,2);
    trans_landmark(:,1) = 641 - trans_landmark(:,1);
end

% Source1.vertices = test_face_neu_new;
% Source1.normals = test_face_neu_new;
end