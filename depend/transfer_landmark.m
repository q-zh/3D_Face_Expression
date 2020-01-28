function land_out = transfer_landmark(land_in)
% input 51*3
% output 5*3

%%
land_out(1, 1:3) = (land_in(20, 1:3) + land_in(23, 1:3))/2;
land_out(2, 1:3) = (land_in(26, 1:3) + land_in(29, 1:3))/2;
land_out(3, 1:3) = land_in(14, 1:3);
land_out(4, 1:3) = land_in(32, 1:3);
land_out(5, 1:3) = land_in(38, 1:3);
% 
% land_out(3:22, :) = land_in(32:end, 1:3);

end