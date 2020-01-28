
function  [vertex, edge] = find_vertex_edge(z, mask, mask_mouth)
% z = flip(z,1);
% mask = flip(mask,1);
% mask_mouth = flip(mask_mouth,1);
%
% z = flip(z,2);
% mask = flip(mask,2);
% mask_mouth = flip(mask_mouth,2);
if nargin == 2
    [gx,gy] = meshgrid(1:size(mask,2), 1:size(mask,1));
    
    gx = gx(mask);
    gy = gy(mask);
    gz = z(mask);
    vx = gx(:);
    vy = gy(:);
    vz = gz(:);
    edge = delaunay(vx,vy);
    
    vertex = [vx,vy, vz];
else
    
    [gx,gy] = meshgrid(1:size(mask,2), 1:size(mask,1));
    mask_nomouth = and(xor(mask_mouth, mask), mask);
    gx1 = gx(mask_nomouth);
    gy1 = gy(mask_nomouth);
    gz1 = z(mask_nomouth);
    gx2 = gx(mask_mouth);
    gy2 = gy(mask_mouth);
    gz2 = z(mask_mouth);
    vx = [gx1(:);gx2(:)];
    vy = [gy1(:);gy2(:)];
    vz = [gz1(:);gz2(:)];
    
    
    edge = delaunay(vx,vy);
    num1 = numel(gx1);
    num2 = numel(gx2);
    for i = 1:num2
        index_v = i+num1;
        for j=1:3
            edge(find(edge(:,j)==index_v),:)=[];
        end
    end
    
    gx = gx(mask_nomouth);
    gy = gy(mask_nomouth);
    gz = z(mask_nomouth);
    vx = gx(:);
    vy = gy(:);
    vz = gz(:);
    
    vertex = [vx,vy, vz];
end



