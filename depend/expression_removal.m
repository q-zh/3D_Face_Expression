function test_face_neu = expression_removal(vertices, label)

%% testing
para_pca_dim_exp = 395;
nose_matrix = repmat(vertices(4659,:), [5094,1]);
vertices = vertices - nose_matrix;

load(['dic_' label '.mat']);

final_out = vertices';
final_out = final_out(:);
DicExp = dictionary(1:para_pca_dim_exp,:);
DicRes = dictionary(para_pca_dim_exp+1:end,:);
norm_dic_exp = sqrt(sum(DicExp.^2));
dic_exp_n = DicExp./repmat(norm_dic_exp,size(DicExp,1),1);
norm_dic = sqrt(sum(dictionary.^2));
dic_n = dictionary./repmat(norm_dic,size(dictionary,1),1);


inputs_face = final_out(:);
inputs = pca_exp_vec' * (inputs_face - mean_exp);
code0 = omp(dic_exp_n'*inputs, dic_exp_n'*dic_exp_n, 60);
test_res_old = DicRes*code0;

for j=1:8
    test = [inputs;test_res_old];
    code = omp(dic_n'*test, dic_n'*dic_n, 60);
    test_res = DicRes*code;
%             outData = DicNeu*code+repmat(meanNeu,1,4);
    disp(['iter error ' num2str(mean(mean(abs(test_res-test_res_old))))]);
    if mean(mean(abs(test_res_old-test_res)))<0.01
        break;
    end
    test_res_old = test_res;
end
test_exp = DicExp*code;
test_face_neu = inputs_face - pca_res_vec*test_res - mean_res;
test_face_neu = reshape(test_face_neu, [3, 5094])';
test_face_neu = test_face_neu + nose_matrix;
end